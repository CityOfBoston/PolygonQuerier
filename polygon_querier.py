# Ethan McIntosh - BPDA Research Division - 2023/02/10
import h3, numpy as np
from matplotlib import pyplot as plt
from matplotlib.path import Path
from matplotlib.patches import PathPatch
from matplotlib.collections import PatchCollection
from shapely import wkt, ops as so, geometry as sg, algorithms as sa
from abc import ABCMeta, abstractmethod
from math import sqrt
from pyproj import Transformer, Proj
from pyproj.aoi import AreaOfInterest
from pyproj.database import query_utm_crs_info

def plot_polygon(ax, poly, **kwargs):
    """Plots a given shapely Polygon on a given matplotlib Axes object. Author credit: Stackoverflow user "mins" 
    (https://stackoverflow.com/questions/55522395/how-do-i-plot-shapely-polygons-and-objects-using-matplotlib),
    who in turn adapted it from geopandas source code 
    (https://github.com/geopandas/geopandas/blob/54fe63d3b9daaca9ce0f50ed5711f78acd56c378/geopandas/plotting.py).
    """
    path = Path.make_compound_path(
        Path(np.asarray(poly.exterior.coords)[:, :2]),
        *[Path(np.asarray(ring.coords)[:, :2]) for ring in poly.interiors])

    patch = PathPatch(path, **kwargs)
    collection = PatchCollection([patch], **kwargs)
    
    ax.add_collection(collection, autolim=True)
    ax.autoscale_view()
    return collection

class PolygonQuerier(metaclass = ABCMeta):
    def __init__(self, wkt_str: str, tolerance_km: float, epsg: int = 4326):
        """
        wkt_str: a Well-known Text (WKT) string representing the latitude/longitude coordinate geometry (WGS84) of the 
        area of interest. This area of interest can either be a polygon or a multipolygon.

        tolerance_km: Permissible level of geographic imprecision in calculating successive poles of inaccessibility, expressed in kilometers.

        epsg (optional): If wkt_str uses a coordinate reference system other than latitude/longitude (EPSG 4326), use this parameter to specify 
        the EPSG code that the wkt_str coordinates use. 
        """
        to_wgs84 = Transformer.from_crs(f"epsg:{epsg}", "epsg:4326", always_xy=True).transform
        self.original_geo = so.transform(to_wgs84, wkt.loads(wkt_str)) # load the WKT geometry in as a Shapely object in WGS84 coordinates
        self.tolerance_km = tolerance_km
        self.total_calls = 0  # keep track of how many calls we make to the user-specified query function
        epsg_meters = self.metric(epsg)
        self.m_to_wgs84 = Transformer.from_crs(f"epsg:{epsg_meters}", "epsg:4326", always_xy=True).transform
        self.wgs84_to_m = Transformer.from_crs("epsg:4326", f"epsg:{epsg_meters}", always_xy=True).transform

    def metric(self, epsg: int) -> int:
        """If the CRS represented by the given EPSG does not use meters as a linear unit, then this method returns the EPSG
        for a UTM zone that intersects the area of interest. If the given EPSG does use meters, this method just returns that EPSG.
        """
        crs_def = Proj(epsg).definition
        i = crs_def.find("units=")
        if i < 0 or crs_def[i + len("units=")].lower() != 'm': # if the specified CRS units are not meters...
            ref = self.original_geo.centroid
            utm_crs_list = query_utm_crs_info(
                datum_name="WGS 84",
                area_of_interest=AreaOfInterest(west_lon_degree=ref.x, south_lat_degree=ref.y, east_lon_degree=ref.x, north_lat_degree=ref.y)
            )
            return utm_crs_list[0].code # ...select the UTM zone that best matches the area of interest to use for metric calculations & plotting
        else: # but if the specified CRS is in meters...
            return epsg # we can use that for metric calculations & plotting

    def run(self, nearest_neighbors_query_func, target_coverage=1.0, max_new_calls=np.Inf, max_total_calls=np.Inf):
        """point_based_query_func: the name of a function that takes in a latitude and longitude and returns a query coverage radius in kilometers.

        target_coverage stops execution of this method once a certain portion of the original area of interest has been covered by query results. Specify any
        value from 0 to 1, with 0 meaning this will perform no queries, and 1 meaning it will perform queries until 100% of the area is covered. If no
        value is specified, the default target_coverage is 1.

        max_new_calls and max_total_calls are used to stop execution based on the number of calls that are made to nearest_neighbors_query_func. Max_new_calls 
        limits how many additional calls are made during a given run() call, while max_total_calls limits the total number of calls made by a given PolygonQuerier.
        Both of these parameters default to infinity (no maximum). 
        """
        new_calls = 0
        if target_coverage > 1: # really this should be flagged as an invalid argument, but for now this is fine
            target_coverage = 1

        while self.coverage() < target_coverage and new_calls < max_new_calls and self.total_calls < max_total_calls:
            self.next_coords_to_query = self.next_coords()
            qlat, qlon = self.next_coords_to_query
            # qcr_km = query coverage radius, in kilometers
            qcr_km = nearest_neighbors_query_func(qlat, qlon)
            # increment our counts of new calls and total calls to point_based_query_func
            new_calls += 1
            self.total_calls += 1

            if qcr_km < 0:
                raise Exception(f"Query coverage radius cannot be a negative number ({qcr_km})")
            if qcr_km == 0:
                print(f"Warning: point {self.next_coords_to_query} returned a query coverage radius of 0. Skipping that point.")
            self.process_covered_radius(qcr_km)

    def coords_in_original_area(self, lat, lon) -> bool:
        """Returns True if the given lat/lon coords fall within the original area of interest, False otherwise."""
        return self.original_geo.contains(sg.Point(lon, lat))
    
    @abstractmethod
    def coords_in_current_area(self, lat, lon) -> bool:
        """Returns True if the given lat/lon coords fall within the original area of interest, False otherwise."""
        pass

    @abstractmethod
    def next_coords(self) -> tuple:
        """Returns a (lat, lon) tuple representing the coordinates of the next point within the area of interest that will be queried."""
        pass

    @abstractmethod
    def coverage(self) -> float:
        """Returns a float between 0 and 1 representing the fraction of the area of interest that has been covered by query results."""
        pass
        
    @abstractmethod
    def process_covered_radius(self, qcr_km: float):
        """Subtracts an approximately circular area of radius qcr_km (query coverage radius in kilometers) from the area of interest centered
        on the most recently queried point, to help keep track of which parts of the area of interest are covered by query results."""
        pass

    @abstractmethod
    def plot(self, save_file : str = None):
        """Plots the area of interest (minus any parts of it that have been covered by query results) onto the console, using projected coordinates. 
        Also plots the point where the next query will be made (the pole of inaccesibility). 

        Optional param save_fig: filepath with extension to save the plot image as a local file
        """
        pass 

class PolylabelPQ(PolygonQuerier):
    def __init__(self, *args, **kwargs):
        super(self.__class__, self).__init__(*args, **kwargs)
        self.current_geo = wkt.loads(self.original_geo.wkt)
        self.next_coords_to_query = self.next_coords()

    def coverage(self) -> float:
        projected_current = so.transform(self.wgs84_to_m, self.current_geo)
        projected_original = so.transform(self.wgs84_to_m, self.original_geo)
        return 1 - projected_current.area / projected_original.area # this currently uses square degrees, which isn't ideal

    def next_coords(self) -> tuple:
        aoi_m = so.transform(self.wgs84_to_m, self.current_geo)

        if aoi_m.geom_type == 'Polygon':
            largest_poly = aoi_m 
        elif aoi_m.geom_type == 'MultiPolygon':
            largest_poly = None
            largest_poly_area = 0
            for geom in aoi_m.geoms:
                if geom.area > largest_poly_area:
                    largest_poly = geom
                    largest_poly_area = geom.area
                    
        return tuple(reversed(so.transform(self.m_to_wgs84, sa.polylabel.polylabel(largest_poly, self.tolerance_km*1000)).coords[0]))

    def process_covered_radius(self, qcr_km: float):
        if qcr_km == 0: # If 0 results come back for a point, we don't want to keep querying it over and over
            buffer_radius_m = self.tolerance_km * 1000 # this feels quick-and-dirty, but I'm not sure how else to handle the 0-result case
        else:
            buffer_radius_m = qcr_km*1000

        aoi_m = so.transform(self.wgs84_to_m, self.current_geo) # project the polygon and the queried point into meters
        point_coords_m = self.wgs84_to_m(self.next_coords_to_query[1],self.next_coords_to_query[0]) 
        buffer_circle = sg.Point(point_coords_m).buffer(buffer_radius_m) # this buffer represents the circular area covered by results
        diff_poly_m = aoi_m.difference(buffer_circle) # subtract the buffer polygon from the area of interest polygon
        self.current_geo = so.transform(self.m_to_wgs84, diff_poly_m) # and convert back to WGS84

    def coords_in_current_area(self, lat, lon) -> bool:
        return self.current_geo.contains(sg.Point(lon, lat))
    
    def plot(self, save_file: str =None):
        fig, ax = plt.subplots()

        if self.current_geo.is_empty:
            print("Area of interest is fully covered by results.")
            return

        # plot the hexagons
        if self.current_geo.geom_type == 'Polygon':
            plot_polygon(ax, so.transform(self.wgs84_to_m, self.current_geo), facecolor='lime')
        elif self.current_geo.geom_type == 'MultiPolygon':
            for geom in self.current_geo.geoms:
                plot_polygon(ax, so.transform(self.wgs84_to_m, geom), facecolor='lime')
        projected_pt = so.transform(self.wgs84_to_m, sg.Point(self.next_coords_to_query[1], self.next_coords_to_query[0]))

        # plot the next point to be queried
        plt.plot(projected_pt.x, projected_pt.y, marker="o", markeredgecolor="red", markerfacecolor="red")

        ax.set_aspect("equal")
        ax.get_xaxis().set_visible(False)
        ax.get_yaxis().set_visible(False)

        if save_file is not None:
            plt.savefig(save_file)

class H3PQ(PolygonQuerier):
    def __init__(self, *args, **kwargs):
        super(self.__class__, self).__init__(*args, **kwargs)

        if self.tolerance_km < h3.edge_length(15, unit='km'):
            raise Exception(f"Specified tolerance is lower than the minimum possible precision of {round(h3.edge_length(15, unit='km'), 6)} km.")
        else:
            for p in range(16):
                if self.tolerance_km >= h3.edge_length(p, unit='km'):
                    self.h3_precision = p
                    break

        self.build_hexes() # construct a set of H3 hex IDs corresponding to the area of interest (self.hexes)
        self.initial_hex_count = len(self.hexes)
        self.compute_initial_hdes() # computes how far away each hex is from the perimeter of the area, assigns self.most_central_hex_id 
        self.next_coords_to_query = self.next_coords()

    def build_hexes(self):
        """Takes the user-input geometry and H3 precision level and initializes self.hexes as a set of H3 cells whose centers are within
        the given area of interest. Raises a TypeError if the provided geometry is not a polygon or a multipolygon.
        """
        if self.original_geo.geom_type == 'Polygon':
            self.hexes = dict.fromkeys(h3.polyfill_geojson(self.original_geo.__geo_interface__, self.h3_precision), None)
        elif self.original_geo.geom_type == 'MultiPolygon':
            hex_ids = set()
            for g in self.original_geo.geoms:
                hex_ids.update(h3.polyfill_geojson(g.__geo_interface__, self.h3_precision))
            self.hexes = dict.fromkeys(hex_ids, None)
        else:
            raise TypeError("WKT geometry must be of type POLYGON(...) or MULTIPOLYGON(...)")

    def compute_initial_hdes(self):
        """HDE is an integer representing a given hex's distance from the edge of the area of interest, counted in numbers of hexes. This method
        calculates HDE for all hexes in the area and returns the hex ID with the highest HDE (the approximate "pole of inacessibility" of the area). 
        
        The algorithm used here for initializing HDEs is simple and functional, but performance degrades quickly as hex precision goes up.
        Luckily, this step needs to be taken only once per PolygonQuerier, and recomputing HDEs after they're initialized (e.g. in response 
        to data source queries) using the recompute_hdes() method is more performant.
        """
        hexes_to_check = set(self.hexes.keys()) # make a copy of the set of hexes representing the area of interest
        hde = 0
        
        while hexes_to_check:
            hexes_to_remove = set()

            for h in hexes_to_check:
                for neighbor in h3.k_ring_distances(h, 1)[1]: # look through each hex's first layer of neighboring hex IDs
                    if neighbor not in hexes_to_check: # if the hex has a neighbor outside the current set of hexes (meaning it's on the "edge")
                        self.hexes[h] = hde # assign it the current iteration's HDE value,
                        hexes_to_remove.add(h) # mark it for removal at the end of this iteration,
                        break # and proceed to the next hex
            
            if len(hexes_to_remove) == len(hexes_to_check): # this is only true when we've reached our maximum HDE
                self.most_central_hex_id = self.choose_most_central(hexes_to_check) # find and assign the "most central" hex among the hexes with the highest HDE
            
            for hex in hexes_to_remove: # remove the latest "edge" of hexes from our copy of the set of hexes
                hexes_to_check.remove(hex)

            hde += 1

    def choose_most_central(self, same_hde_hexes: set) -> str:
        """HDE is an integer representing a given hex's distance from the edge of the area of interest, counted in numbers of hexes.

        Given a set of hex IDs with the same HDE, choose and return the "most central" hex within that set.
        We define that here as the hex with the most neighbors of the same HDE. There may be multiple, in which
        case, we just pick one of those.
        """
        most_central_hex_id = None
        max_neighbor_count = -1

        for hx in same_hde_hexes:
            neighbor_hdes = [self.hexes[h] for h in h3.k_ring_distances(hx, 1)[1] if h in same_hde_hexes]
            same_hde_neighbor_count =  len(neighbor_hdes)

            if same_hde_neighbor_count > max_neighbor_count:
                max_neighbor_count = same_hde_neighbor_count
                most_central_hex_id = hx
        
        return most_central_hex_id

    def recompute_hdes(self, qcr_hx: int):
        """HDE is an integer representing a given hex's distance from the edge of the area of interest, counted in numbers of hexes.
        
        This method makes a single pass through self.hexes, recomputing all of their HDEs and re-assigning self.most_central_hex_id
        based on each hex's distance from the most recently queried point (hd_qp) and the query's coverage radius, expressed in an 
        integer number of hexes (qcr_hx). For example, if a hex is 4 hexes away from the query point and query covered 3 hexes, it is 
        assigned a new HDE value of 0, since it's effectively "on the edge" of the shape formed by cutting out the query coverage circle.
        """
        max_hde_hexes = set()
        max_hde = 0

        for hex_id, hde in self.hexes.items():
            # calculating each hex's distance from the most recently queried point is a simple H3 lookup operation
            hd_qp = h3.h3_distance(self.most_central_hex_id, hex_id)
            new_hde = hd_qp - qcr_hx - 1

            if new_hde >= 0 and new_hde < hde:
                self.hexes[hex_id] = new_hde

            # Meanwhile, we want to keep track of both the new maximum HDE and the hex IDs which have that maximum HDE...
            if self.hexes[hex_id] > max_hde:
                max_hde = self.hexes[hex_id]
                max_hde_hexes = set()
                max_hde_hexes.add(hex_id)
            elif self.hexes[hex_id] == max_hde:
                max_hde_hexes.add(hex_id)

        # ...so that we can re-identify the "most central" hex (i.e. the next point within the area of interest we should query)
        self.most_central_hex_id = self.choose_most_central(max_hde_hexes)

    def coverage(self) -> float:
        return 1 - len(self.hexes) / self.initial_hex_count

    def next_coords(self) -> tuple:
        return h3.h3_to_geo(self.most_central_hex_id)

    def process_covered_radius(self, qcr_km: float):
        # icr_km = inscribed circle radius (in km) of the average hex at the chosen level of precision (geometrically related to average edge length) 
        icr_km = h3.edge_length(self.h3_precision, unit='km') * sqrt(3) / 2 # we don't need to recalculate this each time, but it's not expensive either
        # qcr_hx = query coverage radius, in number of hexes. This has a geometric relationship with qcr_km and icr_km
        qcr_hx = int((qcr_km/icr_km - 1)/2)

        for hxs in h3.k_ring_distances(self.most_central_hex_id, qcr_hx):
            for hx in hxs: # any hexes within the query coverage radius should be removed
                self.hexes.pop(hx, None) # the None argument ensures that no KeyError is raised if hx is not a member of hexes

        self.recompute_hdes(qcr_hx) # recomputes self.most_central_hex_id

    def coords_in_current_area(self, lat, lon) -> bool:
        return h3.geo_to_h3(lat, lon, self.h3_precision) in self.hexes
    
    def plot(self, save_file: str =None):
        if not self.hexes:
            print("Area of interest is fully covered by results.")
            return
        
        # convert the set of hexes into a list of Shapely polygons representing each hex, then calculate the unary union of them all
        hex_polys = []
        for hx in self.hexes.keys():
            xy_coords = []
            for elem in h3.h3_to_geo_boundary(hx): # h3 spits out hex coordinates using (lat, lon) tuples
                xy_coords.append(self.wgs84_to_m(elem[1], elem[0])) # but shapely polygons take (x, y) tuples
            hex_polys.append(sg.Polygon(xy_coords)) # so we flip each tuple when building the shapely polygons
        hxs = so.unary_union(hex_polys) # that way, the resulting object is plotted with lon on the x-axis and lat on the y-axis

        # plot the combined polygon/multipolygon onto the matplotlib Axes object
        fig, ax = plt.subplots()
        if hxs.geom_type == 'Polygon':
            plot_polygon(ax, hxs, facecolor='lime')
        elif hxs.geom_type == 'MultiPolygon':
            for poly in hxs.geoms:
                plot_polygon(ax, poly, facecolor='lime')

        # plot the most central hex as a red point
        y, x = h3.h3_to_geo(self.most_central_hex_id)
        x, y = self.wgs84_to_m(x, y)
        plt.plot(x, y, marker="o", markeredgecolor="red", markerfacecolor="red")

        ax.set_aspect("equal")
        ax.get_xaxis().set_visible(False)
        ax.get_yaxis().set_visible(False)

        if save_file is not None:
            plt.savefig(save_file)