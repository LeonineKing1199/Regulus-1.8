# Regulus-1.8
Degenerate tetrahedralization


Alright, this is something I'm definitely far more proud of than v1.5 which is now horribly outdated.

This code is a modification of the gFlip3D/gDel3D algorithm found at http://www.comp.nus.edu.sg/~tants/gdel3d.html

I think the authors are absolutely amazing so I decided to build on their work and approach the triangulation process a bit differently.

For those not in the know, a very fun thing to do in mathematics is take a set of points and draw triangles between them, creating a mesh. It's become common to approach these problems by creating an all-encompassing triangle (or tetrahedron in 3D) and then inserting points one at a time. Each point insertion creates a "fracture" or a set of sub-tetrahedra.gFlip3D struck me as unique because unlike CPUs which insert points one at a time, the authors used the GPU to insert many points at a time.

However, it's difficult to handle inserting points when they are on a tetrahedron (either on a face or an edge). In fact, it complicates things quite a bit.

A popular way to address these types of points is to "perturb" them. For complete info, see the paper on the Simulation of Simplicity which explains the technique in great detail. I, however, choose to modify the point insertion algorithm itself in lieu of modifying the points and how they are represented. What results is most likely slower code but code that I feel is far more elegant.

Paper coming soon!
