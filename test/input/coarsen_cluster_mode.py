#!/usr/bin/env python3
"""Area-weighted dominant-cluster aggregation from 0.5° to 1.25° x 1.875°.

Usage:
    python coarsen_cluster_mode.py INPUT.nc --output Cluster_1p25x1p875.nc --figure Cluster_comparison.png

Dependencies: numpy, scipy, h5py, matplotlib.
The input must be a NetCDF4/HDF5 file containing 2-D Cluster(lat,lon)
and 1-D lat and lon at 0.5-degree cell centers. Missing cells are ignored.
When two types tie, the smaller numerical cluster ID is selected.
"""
import argparse
import numpy as np
import h5py
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import BoundaryNorm, ListedColormap
from scipy.sparse import csr_matrix
from scipy.io import netcdf_file

FILL = -99999


def overlap_matrix(source_edges, target_edges, latitude=False):
    """Sparse target-by-source fractional intersection weights, with spherical area weighting for latitude."""
    rows, cols, weights = [], [], []
    for j in range(len(target_edges)-1):
        lo, hi = target_edges[j:j+2]
        first = max(0, int(np.searchsorted(source_edges, lo, side='right')-1))
        last = min(len(source_edges)-1, int(np.searchsorted(source_edges, hi, side='left')+1))
        for i in range(first, last):
            a, b = max(lo, source_edges[i]), min(hi, source_edges[i+1])
            if b <= a:
                continue
            w = (np.sin(np.deg2rad(b))-np.sin(np.deg2rad(a))) if latitude else (b-a)
            rows.append(j); cols.append(i); weights.append(w)
    return csr_matrix((weights, (rows,cols)), shape=(len(target_edges)-1,len(source_edges)-1))


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('input', help='Original valid NetCDF4 (.nc) file')
    ap.add_argument('--output', default='Cluster_1p25x1p875.nc')
    ap.add_argument('--figure', default='Cluster_comparison.png')
    args = ap.parse_args()

    with h5py.File(args.input, 'r') as fin:
        clusters = np.asarray(fin['Cluster'][:])
        lat = np.asarray(fin['lat'][:], dtype=float)
        lon = np.asarray(fin['lon'][:], dtype=float)
        fill = int(np.asarray(fin['Cluster'].attrs.get('_FillValue', [FILL])).flat[0])
    if clusters.ndim != 2 or clusters.shape != (lat.size, lon.size):
        raise ValueError('Expected Cluster(lat,lon) 2-D field.')
    if np.any(np.diff(lat) < 0):
        lat = lat[::-1]; clusters = clusters[::-1, :]
    if np.any(np.diff(lon) < 0):
        lon = lon[::-1]; clusters = clusters[:, ::-1]
    if not (np.allclose(np.diff(lat), .5) and np.allclose(np.diff(lon), .5)):
        raise ValueError('Expected regularly spaced 0.5-degree source coordinates.')
    if not (np.isclose(lat[0], -89.75) and np.isclose(lat[-1],89.75)
            and np.isclose(lon[0],-179.75) and np.isclose(lon[-1],179.75)):
        raise ValueError('Expected global coordinates centered from -89.75 to 89.75 and -179.75 to 179.75.')

    lat_edges = np.arange(-90, 90.0001, .5)
    lon_edges = np.arange(-180, 180.0001, .5)
    target_lat_edges = np.linspace(-90,90,145)       # 144 bands of 1.25 degrees
    target_lon_edges = np.linspace(-180,180,193)     # 192 bands of 1.875 degrees
    target_lat = (target_lat_edges[:-1]+target_lat_edges[1:])/2
    target_lon = (target_lon_edges[:-1]+target_lon_edges[1:])/2
    lat_weights = overlap_matrix(lat_edges, target_lat_edges, latitude=True)
    lon_weights = overlap_matrix(lon_edges, target_lon_edges)

    valid = np.isfinite(clusters) & (clusters != fill)
    ids = np.unique(clusters[valid].astype(int))
    if ids.size == 0:
        raise ValueError('No valid Cluster values found.')
    # Compare spherical surface-area contribution of each cluster within each target cell.
    # This correctly accounts for fractional source cells along misaligned boundaries.
    score = np.empty((len(ids),len(target_lat),len(target_lon)), dtype=np.float64)
    for k, cid in enumerate(ids):
        indicator = (valid & (clusters == cid)).astype(float)
        score[k] = (lat_weights @ indicator) @ lon_weights.T
    winners = np.argmax(score, axis=0)  # sorted IDs => smallest ID wins a tie
    maxscore = np.max(score, axis=0)
    coarse = np.where(maxscore > 0, ids[winners], FILL).astype(np.int32)

    with netcdf_file(args.output,'w') as out:
        out.createDimension('lat', len(target_lat))
        out.createDimension('lon', len(target_lon))
        yy = out.createVariable('lat','f4',('lat',))
        xx = out.createVariable('lon','f4',('lon',))
        cc = out.createVariable('Cluster','i4',('lat','lon'))
        yy[:] = target_lat.astype('f4'); xx[:] = target_lon.astype('f4')
        cc[:] = coarse
        yy.units = 'degrees_north'; xx.units = 'degrees_east'
        cc.long_name = 'Dominant K-means Cluster (area-weighted mode)'
        cc._FillValue = np.int32(FILL)
        out.description = 'Global 1.25 degree latitude x 1.875 degree longitude; spherical area-weighted modal cluster from 0.5 degree source. Missing cells ignored; lower ID wins ties.'

    # Use common fixed categorical colors, with missing cells light gray.
    import matplotlib.colors as mcolors
    colors = list(plt.get_cmap('tab10').colors) if len(ids) <= 10 else list(plt.get_cmap('tab20',len(ids)).colors)
    cmap = ListedColormap(colors[:len(ids)]); cmap.set_bad('#e5e5e5')
    boundaries = np.arange(len(ids)+1)-.5
    norm = BoundaryNorm(boundaries,len(ids))
    def category_image(data):
        indices = np.searchsorted(ids, data)
        return np.ma.array(indices, mask=(data == FILL) | ~np.isin(data,ids))
    fig, axs = plt.subplots(1,2,figsize=(16,6.4), constrained_layout=True)
    for ax, data, title in [(axs[0],clusters,'Original: 0.5° × 0.5° (360 × 720)'),
                             (axs[1],coarse,'Coarsened: 1.25° × 1.875° (144 × 192)')]:
        im=ax.imshow(category_image(data),origin='lower',extent=(-180,180,-90,90),
                     interpolation='nearest',aspect='auto',cmap=cmap,norm=norm)
        ax.set(xlabel='Longitude (°)',ylabel='Latitude (°)',title=title,
               xticks=[-180,-120,-60,0,60,120,180],yticks=[-90,-60,-30,0,30,60,90])
        ax.grid(color='black',alpha=.12,linewidth=.5)
    cb=fig.colorbar(im,ax=axs,orientation='horizontal',shrink=.68,pad=.10,
                    ticks=np.arange(len(ids)))
    cb.ax.set_xticklabels([str(c) for c in ids]);cb.set_label('Cluster ID (grey = missing)')
    fig.suptitle('Global soil cluster classification: before and after area-weighted modal aggregation',fontsize=15)
    fig.savefig(args.figure,dpi=220,bbox_inches='tight',facecolor='white')
    plt.close(fig)
    print(f'Input: {clusters.shape}; output: {coarse.shape}; classes: {ids.tolist()}')
    print(f'Wrote {args.output} and {args.figure}')

if __name__ == '__main__':
    main()
