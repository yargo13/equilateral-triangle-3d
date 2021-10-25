import numpy as np

def latLonToSpherical(latLonDegrees):
    sphericalDegrees = np.array([
        90 - latLonDegrees[0],
        180 + latLonDegrees[1]
    ])
    return np.radians(sphericalDegrees)

def sphericalToLatLon(sphericalRadians):
    sphericalDegrees = np.degrees(sphericalRadians)
    return np.array([
        90 - sphericalDegrees[1],
        (sphericalDegrees[2] - 180)%360
    ])

def toSpherical(xyz):
    xy = xyz[0]**2 + xyz[1]**2
    r = np.sqrt(xy + xyz[2]**2)
    theta = np.arctan2(np.sqrt(xy), xyz[2]) 
    phi = np.arctan2(xyz[1], xyz[0])
    return np.array([r, theta, phi])

def toCartesian(rthetaphi):
    return np.array([
        rthetaphi[0]*np.sin(rthetaphi[1])*np.cos(rthetaphi[2]),
        rthetaphi[0]*np.sin(rthetaphi[1])*np.sin(rthetaphi[2]),
        rthetaphi[0]*np.cos(rthetaphi[1])
    ])

def RU(source, destination):
    cross = np.cross(source, destination)
    ssc = np.array([
        [0, -cross[2], cross[1]],
        [cross[2], 0, -cross[0]],
        [-cross[1], cross[0], 0]
    ])
    return np.eye(3) + ssc + \
        ssc.dot(ssc).dot((1-np.dot(source, destination))/(np.linalg.norm(cross)**2))

def rotZ(angle):
    c = np.cos(angle)
    s = np.sin(angle)
    return np.array([
        [c, -s, 0],
        [s, c, 0],
        [0, 0, 1]
    ])


def calculateRotation(source1, source2, destination1):
    RU1 = RU(source1, destination1)
    RUSourceToDestination = np.linalg.inv(RU1)
    possibleDestination = source2.dot(RUSourceToDestination)
    possibleDestinationSpherical = toSpherical(possibleDestination)
    ROTZ = rotZ(possibleDestinationSpherical[2])
    return np.linalg.inv(RUSourceToDestination.dot(ROTZ))

# Coordinates of Point 1
lat1 = 49.2804609178205
lon1 = -123.11971971722261
# Coordinates of Point 2
lat2 = -23.558764358499484
lon2 = -46.68946125518872

latLon1 = np.array([lat1, lon1])
latLon2 = np.array([lat2, lon2])
latLonRadians1 = latLonToSpherical(latLon1)
latLonRadians2 = latLonToSpherical(latLon2)
cartesian1 = toCartesian(np.hstack((1, latLonRadians1)))
cartesian2 = toCartesian(np.hstack((1, latLonRadians2)))

sphericalRotated1 = np.array([1, 0, 0])
alpha = np.arccos(np.dot(cartesian1, cartesian2))
sphericalRotated2 = np.array([1, alpha, 0])
beta = np.arccos(np.cos(alpha)/(np.cos(alpha)+1))
sphericalRotatedResult = np.array([1, alpha, beta])

cartesianRotated1 = toCartesian(sphericalRotated1)
cartesianRotated2 = toCartesian(sphericalRotated2)
cartesianRotatedResult = toCartesian(sphericalRotatedResult)

rotation = calculateRotation(cartesian1, cartesian2, cartesianRotated1)
cartesianResult = cartesianRotatedResult.dot(rotation)
sphericalResult = toSpherical(cartesianResult)
latLonResult = sphericalToLatLon(sphericalResult)
print(latLonResult)
