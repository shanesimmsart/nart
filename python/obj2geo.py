# Python script to convert .obj files to .geo files
# Note that nart uses Z as up and Y as forward
# (Right-handed)

import sys

f = open(sys.argv[1])
geoFile = sys.argv[1].split("obj")[0] + "geo"
g = open(sys.argv[1].split("obj")[0] + "geo", "w+")

fString = f.read()
gString = str(fString.count("f ")) + " \n"

# faces
faceString = ""
normIndexString = ""

n = 0
for item in fString.split("f ")[1:]:
	gString += str(item.split("\n")[0].count(" ")+1) + " "

	for face in item.split("\n")[0].split(" "):
		faceString += str(int(face.split("/")[0])-1) + " "
		normIndexString += str(int(face.split("/")[2])-1) + " "
	n += 1

print("Faces: " + str(n))

gString += "\n" + faceString + "\n"

# verts
vertString = ""

n = 0
for vert in fString.split("v ")[1:]:
	for coord in vert.split("vn")[0].split("vt")[0].split(" "):
		vertString += coord.split("\n")[0] + " "
		# print coord # .split("\n")[0]
	n += 1

print("Vertices: " + str(n))

gString += vertString + "\n"

# normals
normString = ""

n = 0
for vert in fString.split("vn ")[1:]:
	for coord in vert.split("usemtl")[0].split(" "):
		normString += coord.split("\n")[0] + " "
	n += 1
print("Vertex normals: " + str(n))

gString += normIndexString + "\n"

gString += normString

g.write(gString)

print("Created " + geoFile)


