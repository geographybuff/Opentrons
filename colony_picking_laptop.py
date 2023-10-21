import math
from matplotlib import pyplot as plt
import string
import socket

def labwareJson(coords, zDimension, depth, diameter, wellVolume, zBottom, displayName):
    keys = ["ordering", "metadata", "schemaVersion", "version", "namespace", "dimensions", "parameters", "wells", "brand", "groups", "cornerOffsetFromSlot"]
    number_of_colonies = len(coords)
    def Ordering(number_of_colonies):
        letters = list(string.ascii_uppercase)
        well_ids = [[]]
        for i in range(number_of_colonies):
            well_ids[0].append(letters[i%26] + str(i//26+1))
        return well_ids
    ordering = Ordering(number_of_colonies)
    metadata = {
        "displayName": displayName,
        "displayVolumeUnits": "mL",
        "displayCategory": "other",
        "tags": []
    }
    schemaVersion = 2
    version = 1
    namespace = "custom_beta"
    dimensions = {
        "xDimension": 127.76,
        "yDimension": 85.47,
        "zDimension": zDimension
    }
    parameters = {
        "format": "irregular",
        "isTiprack": False,
        "isMagneticModuleCompatible": False,
        "loadName": displayName
    }
    def Wells(coords, depth, diameter, wellVolume, zBottom):
        well_data = {}
        for i in range(number_of_colonies):
            well_data[ordering[0][i]] = {
                "shape": "circular",
                "depth": depth,
                "diameter": diameter,
                "totalLiquidVolume": wellVolume,
                "x": coords[i][0],
                "y": coords[i][1],
                "z": zBottom
            }
        return well_data
    wells = Wells(coords, depth, diameter, wellVolume, zBottom)
    brand = {
        "brand": "Open Insulin",
        "brandId": [],
        "links": ["openinsulin.org", "chat.openinsulin.org", "gitlab.com/open-insulin", "jitsi.openinsulin.org/openinsulin"]
    }
    groups = [{
        "wells": ordering[0],
        "metadata": {
            "displayCategory": "other",
            "wellBottomShape": "flat"
        }
    }]
    cornerOffsetFromSlot = {
        "x": 0,
        "y": 0,
        "z": 0
    }
    values = [ordering, metadata, schemaVersion, version, namespace, dimensions, parameters, wells, brand, groups, cornerOffsetFromSlot]
    return str(dict(zip(keys, values))).replace("False", "false").replace("'", '"')

my_results = "e_coli_5_5.csv"

colonies = pd.read_csv(my_results)
colonies = colonies[(colonies["Y"]<2430) & (colonies["X"]<3300)]
colonies = colonies[(colonies["Y"]>500) & (colonies["X"]>500)]

x = list(colonies["X"])
y = list(colonies["Y"])

proximity_threshold = 1
# Proximity of colonies that will be tolerated. 1 means colonies will be rejected
# if they are touching each other. <1 means some colonies touching each other
# will be accepted. >1 means colonies will be rejected if they are too close without touching.
radius = list(colonies["Radius"])
nearest_neighbor = []
nn_index = []
upper_left_sum = x[0]+y[0]  # used to find which colony is furthest to the upper left
upper_left_index = 0
lower_right_sum = x[0]+y[0]  # used to find which colony is furthest to the bottom right
lower_right_index = 0
coords = []
for i in range(len(x)):
    i_neighbors = []
    for j in range(len(x)):
        if i != j:
            i_neighbors.append(math.sqrt((x[i]-x[j])**2+(y[i]-y[j])**2))
    nearest = min(i_neighbors)
    nearest_neighbor.append(nearest)
    nn_index.append(i_neighbors.index(nearest))
    if radius[i]>=14 and nearest_neighbor[i]>(radius[i]+radius[nn_index[i]])*proximity_threshold:
        coords.append((x[i], y[i]))
        if x[i]+y[i]<upper_left_sum and radius[i]>=14:
            upper_left_sum = x[i]+y[i]
            upper_left_index = i
        if x[i]+y[i]>lower_right_sum and radius[i]>=14:
            lower_right_sum = x[i]+y[i]
            lower_right_index = i
# the indices of the extreme colonies go in the last item when sending to OT-2 
filtered_upper_left_index = coords.index((x[upper_left_index], y[upper_left_index]))
filtered_lower_right_index = coords.index((x[lower_right_index], y[lower_right_index]))
coords.append((filtered_upper_left_index, filtered_lower_right_index))
extreme_indices = [upper_left_index, lower_right_index]
x_without_extremes = x[:min(extreme_indices)]+x[min(extreme_indices)+1:max(extreme_indices)]+x[max(extreme_indices)+1:]
y_without_extremes = y[:min(extreme_indices)]+y[min(extreme_indices)+1:max(extreme_indices)]+y[max(extreme_indices)+1:]
plt.scatter(x_without_extremes, y_without_extremes, c="b")
plt.scatter([x[upper_left_index]], [y[upper_left_index]], c="r")
plt.scatter([x[lower_right_index]], [y[lower_right_index]], c="lime")
plt.xlabel("X coordinate")
plt.ylabel("Y coordinate")
plt.title("Red: Furthest upper left colony; Green: Furthest lower right colony")
plt.ylim(max(y)+80, min(y)-80)
plt.show()

ulname = "ul"+str(x[upper_left_index])+"x"+str(y[upper_left_index])+"y.json"
lrname = "lr"+str(x[lower_right_index])+"x"+str(y[lower_right_index])+"y.json"
ulLabware = labwareJson([(0,0)], zDimension = 25, depth = 5, diameter=5, wellVolume = 200, zBottom = 20, displayName = ulname[:-5])
lrLabware = labwareJson([(0,0)], zDimension = 25, depth = 5, diameter=5, wellVolume = 200, zBottom = 20, displayName = lrname[:-5])
labware_folder = "C:/Users/joshu/AppData/Roaming/Opentrons/labware/"
ulfile = open(labware_folder + ulname, "w")
ulfile.write(ulLabware)
ulfile.close()
lrfile = open(labware_folder + lrname, "w")
lrfile.write(lrLabware)
lrfile.close()

clslot = int(input("Enter slot number for colony plate: "))
trslot = int(input("Enter tiprack slot: "))


calibration_template = open("calibration_template.txt", "r")
calibration_content = calibration_template.read()
calibration_content = calibration_content.replace("cl_name", ulname[:-5]).replace("cl_slot", str(clslot))
calibration_content = calibration_content.replace("tiprack_slot", str(trslot))
calibration_filename = "calibration"+ulname[:-5]+".py"
calibration_file = open("calibration"+ulname[:-5]+".py", "w")
calibration_file.write(calibration_content)
calibration_file.close()
print("Run the calibration file entitled " + calibration_filename + ". Click on Run Labware Position Check and calibrate to the upper left colony.")
ulx = float(input("Enter x calibration offset for " + ulname[:-5] + ": "))
uly = float(input("Enter y calibration offset for " + ulname[:-5] + ": "))
ulz = float(input("Enter z calibration offset for " + ulname[:-5] + ": "))
calibration_template = open("calibration_template.txt", "r")
calibration_content = calibration_template.read()
calibration_content = calibration_content.replace("cl_name", lrname[:-5]).replace("cl_slot", str(clslot))
calibration_content = calibration_content.replace("tiprack_slot", str(trslot))
calibration_filename = "calibration"+lrname[:-5]+".py"
calibration_file = open("calibration"+lrname[:-5]+".py", "w")
calibration_file.write(calibration_content)
calibration_file.close()
print("Run the calibration file entitled " + calibration_filename + ". Click on Run Labware Position Check and calibrate to the upper left colony.")
lrx = float(input("Enter x calibration offset for " + lrname[:-5] + ": "))
lry = float(input("Enter y calibration offset for " + lrname[:-5] + ": "))
lrz = float(input("Enter z calibration offset for " + lrname[:-5] + ": "))

message = coords + [(ulx, lrx), (uly, lry), (ulz, lrz), (clslot, clslot), (trslot, 0)]

HOST = "169.254.40.170"  # The server's hostname or IP address
PORT = 48889  # The port used by the server

with socket.socket(socket.AF_INET, socket.SOCK_STREAM) as s:
    s.connect((HOST, PORT))
    #if len()
    s.sendall(str(message).encode("ASCII"))
    #coordinates must be encoded before sending to ensure cross compatibility
    data = s.recv(8192)