import folium
import pandas as pd
import geopandas as gpd
import pyproj
from shapely.geometry import Point, Polygon
from shapely.ops import transform
from folium.plugins import MeasureControl, Draw, MarkerCluster

# Geometry creation functions
def create_geodesic_shape(center, radius, angle=None, beamwidth=120):
    geod = pyproj.Geod(ellps='WGS84')
    points = []
    
    if angle is None:  # Create a circle
        for i in range(361):
            lon, lat, _ = geod.fwd(center.x, center.y, i, radius)
            points.append(Point(lon, lat))
    else:  # Create a wedge
        for i in range(int(beamwidth) + 1):
            theta = angle - (beamwidth / 2) + i
            lon, lat, _ = geod.fwd(center.x, center.y, theta, radius)
            points.append(Point(lon, lat))
        # Close the polygon with geodesic lines back to the center
        points.extend([Point(*geod.fwd(center.x, center.y, angle + (beamwidth / 2), radius)[:2]),
                       center,
                       Point(*geod.fwd(center.x, center.y, angle - (beamwidth / 2), radius)[:2])])
    
    return Polygon(points)

# Projection functions
def to_projected(geom, src_crs, dst_crs):
    project = pyproj.Transformer.from_crs(src_crs, dst_crs, always_xy=True).transform
    return transform(project, geom)

def get_utm_crs(center_lon):
    return pyproj.CRS.from_proj4(
        f"+proj=utm +zone={int((center_lon + 180) / 6) + 1} +datum=WGS84 +units=m +no_defs"
    )

# Data processing functions
def load_and_process_data(file_path):
    df = pd.read_csv(file_path)
    gdf = gpd.GeoDataFrame(df, geometry=gpd.points_from_xy(df.longitude, df.latitude))
    gdf['beamwidth'] = gdf['beamwidth'].fillna(120)
    gdf['coverage_area'] = gdf.apply(lambda row: 
        create_geodesic_shape(row.geometry, row['range'], row['sector_id'], row['beamwidth']) 
        if pd.notnull(row['sector_id']) 
        else create_geodesic_shape(row.geometry, row['range']), 
        axis=1)
    return gdf

def project_geometries(gdf, utm_crs):
    gdf['projected_geometry'] = gdf.geometry.apply(lambda geom: to_projected(geom, "EPSG:4326", utm_crs))
    gdf['projected_coverage_area'] = gdf.coverage_area.apply(lambda geom: to_projected(geom, "EPSG:4326", utm_crs))
    return gdf

# Wedge grouping function
def group_coverage_areas(coverage_areas):
    groups = []
    while coverage_areas:
        base = coverage_areas.pop(0)
        group = [base]
        overlaps = [area for area in coverage_areas if base.intersects(area)]
        while overlaps:
            overlap = overlaps.pop(0)
            group.append(overlap)
            coverage_areas = [area for area in coverage_areas if area != overlap]
            overlaps.extend([area for area in coverage_areas if overlap.intersects(area) and area not in overlaps])
        groups.append(group)
    return groups

# Map creation functions
def create_base_map(gdf):
    bounds = gdf.total_bounds
    center_lat, center_lon = (bounds[1] + bounds[3]) / 2, (bounds[0] + bounds[2]) / 2
    
    m = folium.Map(location=[center_lat, center_lon], zoom_start=None)
    m.fit_bounds([[bounds[1], bounds[0]], [bounds[3], bounds[2]]])
    
    m.add_child(MeasureControl())
    m.add_child(Draw(
        draw_options={'polyline': True, 'rectangle': True, 'polygon': True, 'circle': True, 'marker': True, 'circlemarker': False},
        edit_options={'edit': True}
    ))
    return m

def add_coverage_areas_and_intersections(m, groups, utm_crs):
    colors = ['blue', 'green', 'yellow', 'purple', 'orange']
    for group_index, group in enumerate(groups):
        intersection = group[0]
        for area in group[1:]:
            intersection = intersection.intersection(area)
        
        intersection = to_projected(intersection, utm_crs, "EPSG:4326")
        
        for i, area in enumerate(group):
            area_wgs84 = to_projected(area, utm_crs, "EPSG:4326")
            shape_type = 'Circle' if area_wgs84.geom_type == 'Polygon' and len(area_wgs84.exterior.coords) > 100 else 'Wedge'
            folium.GeoJson(
                data=area_wgs84.__geo_interface__,
                style_function=lambda _, color=colors[(group_index * len(group) + i) % len(colors)]: {'fillColor': color, 'color': color, 'fillOpacity': 0.3},
                tooltip=f"{shape_type} {group_index * len(group) + i + 1}"
            ).add_to(m)
        
        folium.GeoJson(
            data=intersection.__geo_interface__,
            style_function=lambda _: {'fillColor': 'red', 'color': 'red', 'fillOpacity': 0.5}
        ).add_to(m)

def add_tower_markers(m, gdf):
    marker_cluster = MarkerCluster().add_to(m)
    for idx, row in gdf.iterrows():
        folium.Marker(location=[row['latitude'], row['longitude']], popup=f"Tower {idx+1}").add_to(marker_cluster)

# Main execution
def main():
    gdf = load_and_process_data('towers.csv')
    utm_crs = get_utm_crs(gdf['longitude'].mean())
    gdf = project_geometries(gdf, utm_crs)
    
    coverage_areas = list(gdf['projected_coverage_area'])
    groups = group_coverage_areas(coverage_areas)
    
    gdf['coverage_area'] = gdf.projected_coverage_area.apply(lambda geom: to_projected(geom, utm_crs, "EPSG:4326"))
    
    m = create_base_map(gdf)
    add_coverage_areas_and_intersections(m, groups, utm_crs)
    add_tower_markers(m, gdf)
    
    m.save('map.html')

if __name__ == "__main__":
    main()