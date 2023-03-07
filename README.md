# PolygonQuerier
#### A Python utility for systematically extracting point-based spatial data within a polygon/multipolygon area of interest from data sources that only support nearest-neighbor queries.

<p align="center">
  <img src="https://github.com/CityOfBoston/PolygonQuerier/blob/main/img/pq_demo.gif">
</p>

Many sources of point-based spatial data do not natively support queries that return all data points within a polygon/multipolygon area of interest. In particular, web-based APIs for point-of-interest location data (e.g. from companies like Google, Bing, Yelp, etc.) tend to support only nearest-neighbor queries, returning a subset of points of interest that are nearest to a given set of coordinates. This works well for applications like navigation, routing, and web search. However, for other applications, it's desirable to pull all available point-of-interest data from within a bounded area of interest, not just for a circular area around a single point.

PolygonQuerier is a tool that executes nearest-neighbor queries iteratively until the circular-area results from those queries sufficiently cover an area of interest. In order to do this efficiently, PolygonQuerier repeatedly calculates the pole of inaccessibility of the area of interest (the point within it that is farthest away from any edge) as it gets covered by results. The above animation shows how PolygonQuerier successively calculates poles of inaccessibility within the city of Boston in order to gradually cover the entire city with results from a point-of-interest data source using nearest-neighbor queries. 

## Setup

PolygonQuerier is not currently published on any Python package managers like PyPI or conda. To use it in your code, download or clone this repository locally. You will also need to have the numpy, matplotlib, pyproj, shapely, and h3 packages installed on your version of Python 3.x for PolygonQuerier to function.

To initialize a PolygonQuerier, you must pass in a WKT string representing your area of interest. You must also specify a tolerance in kilometers indicating the maximum acceptable degree of imprecision when calculating successive poles of inaccessibility. An optional third parameter allows you to specify the integer EPSG code of the coordinate reference system that your WKT string uses. By default, PolygonQuerier assumes that your WKT string uses latitude/longitude (EPSG:4326) coordinates.
```
from polygonquerier.polygon_querier import *

my_area_wkt_string = ... # PolygonQuerier supports both POLYGON and MULTIPOLYGON geometries
my_tolerance = 0.1 # this is in kilometers
my_epsg = 26986
```

You have a choice between two flavors of PolygonQuerier, which use different algorithms for calculating poles of inaccessibility. 

### PolylabelPQ

PolylabelPQ is built on the [polylabel algorithm](https://github.com/mapbox/polylabel) by Volodymyr Agafonkin. This algorithm is precise and should be sufficiently performant for most uses of PolygonQuerier.
```
pq = PolylabelPQ(my_area_wkt_string, my_tolerance)
```

### H3PQ

H3PQ offers an alternative and more approximate algorithm for computing poles of inaccessibility. By doing some extra up-front calculations to represent your area of interest as a set of hexagonal grid cells from the the [H3 spatial indexing system](https://h3geo.org/), H3PQ can reduce total processing times for applications where large numbers of queries are required in order to fully cover the area of interest. For more detail on the algorithms behind H3PQ, read the "Notes on pole-of-inaccessibility algorithms" section of this document.
```
pq = H3PQ(my_area_wkt_string, my_tolerance)
```
Note that when using H3PQ, the effective precision of pole-of-inaccessibility calculations may not exactly match your specified tolerance, owing to the fact that the H3 system uses 16 discrete [levels of precision](https://h3geo.org/docs/core-library/restable/). H3PQ will choose a level of precision that's guaranteed to be equal or better than your specified tolerance. For example, an H3PQ with a specified tolerance of 0.1 km will use hexagons at H3 resolution 10, which have an average edge length globally of around 0.07 km. However, there is a minimum effective precision for H3PQ, owing to the fact that the highest H3 resolution is 15, which has an average edge length of around 0.0005 km (0.5 meters). Specifying a tolerance smaller than this when initializing an H3PQ will raise an exception.

## Running a PolygonQuerier

PolygonQueriers do not handle any of the details of querying particular data sources. Therefore, to run a PolygonQuerier on a particular data source, you must write a function that...
- takes in a latitude and a longitude as arguments
- queries the data source for points of interest near the given latitude & longitude, sorted by distance, such that the returned results cover a circular area
- returns the distance (in kilomters) of the radius of the circle that covers the results that were returned. 
... and pass in that function as an argument to the PolygonQuerier run() method, as demonstrated below:
```
def your_query_function(latitude, longitude, *your custom arguments*):
    ...
    ... # define how you'll query the data source, handle the results that are returned, calculate the coverage radius, etc.
    ...
    
    return query_coverage_radius_km
    
pq.run(lambda lat, lon: your_query_function(lat, lon, *your custom arguments*))
```
pq.run() expects a function with exactly two arguments (lat and lon), so using a lambda in this way allows your_query_function() to be written with other arguments that may be desirable for dealing with your particular data source.

If at any point your_query_function() returns a query coverage radius equal to 0, a warning will be printed in the console. Additionally, a small area approximately equal to a circle centered on that point with your specified tolerance as radius will be considered "covered by results" so that the PolygonQuerier can move along without getting caught in an infinite loop. If your_query_function() ever returns a query coverage radius less than 0, the PolygonQuerier will raise an Exception. 

### Specifying a target coverage

By default, the run() method will call your_query_function() until 100% of the area of interest has been covered by results. However, it may be desirable to only run the PolygonQuerier until a certain portion of the area of interest is covered, since nearest-neighbor queries on ever-smaller sections of the area of interest will yield diminishing returns in terms of new results. Specifying a target_coverage fraction as an argument of run() will make the PolygonQuerier stop running once that fraction of the original area of interest is covered by results:
```
# run your_query_function() until 95% of the area of interest is covered by results
pq.run(
    lambda lat, lon: your_query_function(lat, lon, *your custom arguments*)
    , target_coverage = 0.95
)
```

Note that PolylabelPQ and H3PQ handle target_coverage in slightly different ways. For PolylabelPQ, after each call of your_query_function(), the query_coverage_radius_km is used to subtract a circular buffer area representing the query coverage area from the area of interest, and the resulting geometry's area is divided by the area of the original geometry and compared to target_coverage in order to determine whether to continue calling your_query_function().

For H3PQ, the coverage fraction is calculated by dividing the number of hexagons not currently covered by query results by the original number of hexagons approximating the area of interest. H3PQ uses a conservative approximation to determine which hexagons are covered by the query: only hexagons whose shortest hexagon-distance from the queried hexagon is less than the query coverage radius are removed, as illustrated in the diagram below:

<p align="center">
  <img src="https://github.com/CityOfBoston/PolygonQuerier/blob/main/img/hex_coverage.png">
</p>

The set of hexagons that is removed from consideration will therefore always be a smaller area than the circular query coverage area. This guarantees that none of the original area of interest will be falsely marked as covered by query results, although it does introduce opportunities for duplicate results to be pulled in the space between the circular coverage area and the removed hexagons.

### Specifying maximum calls

Another way of controlling the execution of pq.run() is to directly limit the number of calls that it makes to your_query_function(). You can specify the max_total_calls that a given PolygonQuerier can make during its lifetime, and alternatively (or simultaneously), you can specify max_new_calls for a single run() call, in order to run the PolygonQuerier in a step-by-step fashion. The example below also uses the PolygonQuerier plot() method to display which parts of the area of interest are covered by results after a single call to your_query_function():
```
pq.run(
    lambda lat, lon: your_query_function(lat, lon, *your custom arguments*)
    , target_coverage = 0.95
    , max_total_calls = 100
    , max_new_calls = 1
)
pq.plot()
```

## Additional functionalities

#### pq.plot()

The PolygonQuerier plot() method has an optional string argument save_file which allows you to locally save a given plot figure, using any image file extension supported by the [matplotlib savefig() method](https://matplotlib.org/stable/api/_as_gen/matplotlib.pyplot.savefig.html):
```
pq.plot(save_file = 'your_image.png')
```

#### pq.coords_in_original_area() and pq.coords_in_current_area()

These PolygonQuerier methods take exactly two arguments, a latitude and a longitude, and return a boolean representing whether those coordinates are within the area of interest. These methods are identical, but the former considers the original area of interest, while the latter considers the "current" area of interest (the space that hasn't already been covered by results by the PolygonQuerier). 

These methods offer the ability to better control how results are handled inside of your_query_function(). For example, you can choose one of these methods to filter your results before using them for something else:
```
def your_query_function(latitude, longitude, *your custom arguments*):
    # ...query the data source and store them in a variable called results...
    filtered_results = [r for r in results if pq.coords_in_current_area(r['latitude'], r['longitude'])] 
    # ...do stuff with filtered_results...
    # ...calculate the coverage radius...
    return query_coverage_radius_km
    
pq.run(lambda lat, lon: your_query_function(lat, lon, *your custom arguments*))
```

## Notes on coordinate systems

Regardless of the coordinate system you specify when initializing a PolygonQuerier, the original_geo attribute will be a Shapely polygon using latitude/longitude (EPSG:4326) coordinates, allowing for easier interfacing with data sources, many of which use latitude and longitude. However, PolygonQuerier calculations that involve distance and area (e.g. for keeping track of which areas are covered by query results) will perform on-the-fly conversions of geometries to a projected coordinate system in order to reduce error. The plot() method will also use a projected coordinate system when displaying geometries in order to reduce distortion. 

If you initialize a PolygonQuerier using latitude & longitude or any other coordinate system that does not use meters as a linear unit, it will identify a UTM zone that is centered around your area of interest and use that UTM zone's coordinates for calculations and plotting. However, if you specify an EPSG that uses meters when setting up the PolygonQuerier, it will use your coordinate system for those functions. 

## Notes on pole-of-inaccessibility algorithms

Details on the algorithm used by PolylabelPQ for pole of inaccessibility calculations are available at https://github.com/mapbox/polylabel

The H3PQ implementation of PolygonQuerier performs approximate pole of inaccessibility calculations by representing the area of interest as a set of hexagonal grid cells from the [H3 spatial indexing system](https://h3geo.org/). The initial set of hexagons is built using H3's polygon_to_cells (polyfill) method, which returns the set of hexagons at a given resolution whose centroids are within the area of interest.

The initial pole of inaccessibility of the area of interest is calculated by computing each hexagon's HDE (Hexagon-Distance from the Edge). A given hexagon's HDE can be thought of as the minimum number of adjacent hexagons one must traverse in order to reach a hexagon on the edge of the area of interest. The compute_initial_hdes() method assigns an HDE of 0 to any hexagon for which at least one of its neighbors is not a member of the set of hexagons representing the area of interest. Then, it assigns an HDE of 1 to any remaining hexagon for which at least one neighbor has an HDE of 0. The method continues recursively in this way until all hexagons have been assigned an HDE, with the approximate pole of inaccesibility being the center of the hexagon with the maximum HDE. Ties among multiple hexagons having the maximum HDE are broken by considering the average HDE of each hexagon's neighbors (see the choose_most_central() method). This algorithm can be slow for relatively large areas of interest with relatively high levels of precision because of the multiple passes it makes through the set of hexagons.

Where H3PQ achieves its high level of performance is in the subsequent recalculation of the pole of inaccessibility in reponse to each query, taking advantage of H3 indexing functions as well as the unique geometric properties of hexagonal grids. Recomputing each hexagon's HDE in response to each call to your_query_function() requires only a comparison between its current HDE, its hexagon-distance from the hexagon that was queried (a simple H3 lookup operation), and the hexagon-distance of the radius of the circle marking the area covered by query results. The recompute_hdes() method can therefore identify the new pole of inaccesibility by making a single pass through the set of not-yet-queried hexagons doing constant-time operations.
