'''
############################################################################################
## BUREAU OF METEOROLOGY                                                                    
## AUSTRALIAN CLIMATE SERVICE
## DROUGHT AND CHANGES IN ARIDITY HAZARD TEAM
##
## DATE:             Nov-2025
## SCRIPT:           plotting.py
## AUTHOR:           jessica.bhardwaj@bom.gov.au
##
## DESCRIPTION:      Script for plotting flash drought decomposition outputs.
##
############################################################################################
'''

def mask_ocean(input_xr):
    """
    Inputs:
    - input_xr: input xarray to be masked
    Returns:
    Masked input array
    """
    import regionmask
    import numpy as np
    land_mask = regionmask.defined_regions.natural_earth_v5_0_0.land_10.mask(input_xr)
    masked_xr = input_xr.where(land_mask==0, np.nan)
    return masked_xr

def create_multiplot(nrows, ncols, plot_fontsize, title_size, figsize, plot_extent, plot_marker, plot_list, shapefile, shapefile_linewidth, col_titles, col_titles_all_rows_switch, cmap_list, cbar_labels, cbar_nticks, cbar_tick_labels, cbar_x_y_loc, cbar_extensions, cbar_width, cbar_height, cbar_mins=None, cbar_maxs=None):
    """
    Inputs:
    - nrows: integer number of rows to plot
    - ncols: integer number of columns to plot
    - plot_fontsize: integer number of default fontsize
    - title_size: integer number of title fontsize
    - figsize: tuple of figure dimensions 
    - plot_extent: list of lat/lon mins and maxs in format [lon_min, lon_max, lat_min, lat_max]
    - plot_marker: tuple of lat/lon (optional)
    - plot_list: list of data to plot in format [[row1_data], [row2_data], etc]
    - shapefile: str file location of shapefile
    - shapefile_linewidth: int linewidth of shapefile plot feature
    - col_titles: list of string labels for each column 
    - col_titles_all_rows_switch: switch for [True] label each col in each row or [False] label each col in the first row only.
    - cmap_list: list of cmaps for each row e.g. ['viridis', 'plasma', etc]
    - cbar_labels: list of string labels for each row's colourbar
    - cbar_nticks: integer number of ticks for cbar
    - cbar_tick_labels: list of labels for cbar ticks (if not numerical)
    - cbar_x_y_loc: tuple for starting xy position for cbar
    - cbar_extensions: list of string extension information for each row's colourbar e.g. ['max', 'both', etc]
    - cbar_width: float width of each cbar
    - cbar_height: float height of each cbar
    - cbar_nticks: integer number of ticks on cbar
    - cbar_mins: list of integer vmins for each row, specify None if you want default e.g. [0, 2, etc]
    - cbar_max: list of integer vmaxs for each row, specify None if you want default e.g [10, 5, etc] 
    Returns:
    Multi plot cartopy plot with one colourbar to the left of each row. 
    """
    import matplotlib.pyplot as plt
    import cartopy.crs as ccrs
    import cartopy.feature as cfeature
    from cartopy.feature import ShapelyFeature
    from matplotlib.ticker import MaxNLocator
    import numpy as np 
    import geopandas as gpd
    
    plt.rcParams.update({'font.size': plot_fontsize})
    fig, axes = plt.subplots(nrows=nrows, ncols=ncols, figsize=figsize, subplot_kw={'projection': ccrs.PlateCarree(central_longitude=180)})
    if isinstance(axes, np.ndarray):
        axes = axes.flatten()
    else:
        axes = np.array([axes])
    
    if shapefile == 'NRM':
        shapefile_feature = ShapelyFeature(gpd.read_file('/g/data/ia39/aus-ref-clim-data-nci/shapefiles/data/nrm_regions/nrm_regions.shp').dissolve(by="ClusterAb").geometry, ccrs.PlateCarree())
    else:
        shapefile_feature = ShapelyFeature(gpd.read_file(shapefile).geometry, ccrs.PlateCarree()) if shapefile != None else None
    
    for row_idx, row_data in enumerate(plot_list):
        for col_idx, data in enumerate(row_data):
            ax = axes[row_idx * ncols + col_idx]
            cmap = cmap_list[row_idx] if cmap_list and isinstance(cmap_list, list) else None 
            im = data.plot(ax=ax, transform=ccrs.PlateCarree(), cmap=cmap, vmin=cbar_mins[row_idx], vmax=cbar_maxs[row_idx], add_colorbar=False)
            if plot_marker != None:
                ax.plot(plot_marker[1], plot_marker[0], marker='*', color='black', markersize=15, transform=ccrs.PlateCarree(), zorder=10)
            ax.add_feature(cfeature.COASTLINE)
            ax.add_feature(cfeature.BORDERS, linestyle=':')
            if shapefile != None:
                ax.add_feature(shapefile_feature, facecolor='none', edgecolor='black', linewidth=shapefile_linewidth, zorder=6)
            ax.set_title(f'{col_titles[col_idx]}' if (col_titles_all_rows_switch or row_idx == 0) else '', fontweight='bold', fontsize=title_size, pad=10)
            ax.set_extent(plot_extent)

        cbar_ax = fig.add_axes([cbar_x_y_loc[0], cbar_x_y_loc[1] + cbar_x_y_loc[2] * (nrows-1-row_idx), cbar_width, cbar_height]) if nrows != 1 else fig.add_axes([cbar_x_y_loc[0], cbar_x_y_loc[1], cbar_width, cbar_height])
        # cbar_ax = fig.add_axes([0.92, 0.05 + 0.45 * (nrows-1-row_idx), cbar_width, 0.4])
        cbar = fig.colorbar(im, cax=cbar_ax, orientation='vertical', label=cbar_labels[row_idx] if cbar_labels else '', extend=cbar_extensions[row_idx], ticks=MaxNLocator(nbins=cbar_nticks[row_idx]))
        if cbar_tick_labels != None:
            cbar.ax.set_yticklabels(cbar_tick_labels)


    plt.tight_layout(pad=1.0)
    plt.subplots_adjust(right=0.9)
    plt.show()
