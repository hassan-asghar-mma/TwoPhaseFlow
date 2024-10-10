// Define a box shaped geometry to export a surface mesh as input for
// cfMesh or snappyHexMesh
SetFactory("OpenCASCADE");

// A box is defined by a corner point (x, y, z) and its extensions
// in each direction (x, y, z)
vbox = newv;
Box(vbox) = {0, 0, 0, 5.0e-3, 5.0e-3, 2.5e-3};

// Name patches
Physical Surface("front") = {3};
Physical Surface("back") = {4};
Physical Surface("left") = {1};
Physical Surface("right") = {2};
Physical Surface("top") = {6};
Physical Surface("bottom") = {5};

// Mesh settings
Mesh.MeshSizeMax = 2.0e-3;
Mesh.StlOneSolidPerSurface = 2;
