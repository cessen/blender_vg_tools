#=========================== BEGIN MIT LICENSE BLOCK ==========================
#
# Copyright (c) 2012 Nathan Vegdahl
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to
# deal in the Software without restriction, including without limitation the
# rights to use, copy, modify, merge, publish, distribute, sublicense, and/or
# sell copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in
# all copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
# FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS
# IN THE SOFTWARE.
#
#============================ END MIT LICENSE BLOCK ===========================

bl_info = {
    "name": "Vertex Group Tools",
    "author": "Nathan Vegdahl",
    "version": (0, 2),
    "blender": (2, 6, 5),
    "api": 50557,
    "location": "",
    "wiki_url": "",
    "tracker_url": "",
    "category": "Rigging"}

import bpy
import bmesh
from mathutils import Vector


class VCoordHash:
    """ A spatial hash of the vertices of a mesh, for fast look-up by coordinates.
    """
    def __init__(self, mesh):
        """ Initializes the hash.  Right now it uses some bogus
            hueristics to determine the resolution of the hash.  But it seems
            to work okay.  It's certainly far-and-beyond faster than the
            brute-force that would be done without this hash.
        """
        self.mesh = mesh

        # Determine the bounding box of the mesh
        self.mins = [1000000000.0, 1000000000.0, 1000000000.0]
        self.maxs = [-1000000000.0, -1000000000.0, -1000000000.0]
        for v in mesh.vertices:
            for i in range(3):
                if v.co[i] > self.maxs[i]:
                    self.maxs[i] = v.co[i]
                if v.co[i] < self.mins[i]:
                    self.mins[i] = v.co[i]

        # Determin thee resolution of the hash
        quant = 3 * len(mesh.vertices) ** (1 / 3.0)  # Tweak this for different hash granularities
        self.muls = [0, 0, 0]
        for i in range(3):
            n = abs(self.maxs[i] - self.mins[i])
            if n != 0.0:
                self.muls[i] = (1.0 / n) * quant
            else:
                self.muls[i] = 0.0
        avg = (self.muls[0] + self.muls[1] + self.muls[2]) / 3
        self.muls = [avg, avg, avg]

        # Create the hash.  Every vertex gets recorded in the cell it's in
        # as well as all neighboring cells.  This means that every vertex
        # is recorded in a total of 18 hash cells, not just one.  This overlap
        # is what allows us to guarantee correct "closest vert" search results.
        self.hash = {}
        for v in mesh.vertices:
            q = self.quant_coord(v.co)
            for x in [q[0] - 1, q[0], q[0] + 1]:
                for y in [q[1] - 1, q[1], q[1] + 1]:
                    for z in [q[2] - 1, q[2], q[2] + 1]:
                        if (x, y, z) not in self.hash:
                            self.hash[(x, y, z)] = [v.index]
                        else:
                            self.hash[(x, y, z)] += [v.index]

    def quant_coord(self, coord):
        """ Quantizes the coordinates for hash-lookup.
        """
        x = round(coord[0] * self.muls[0])
        y = round(coord[1] * self.muls[1])
        z = round(coord[2] * self.muls[2])
        return (x, y, z)

    def get_verts(self, coord):
        """ Returns a list of vertex indices in the viscinity
            of the given coords.
        """
        q = self.quant_coord(coord)
        vs = []
        if q in self.hash:
            vs += self.hash[q]
        return vs

    def closest_vert(self, coord):
        """ Returns the index of the closest vertex
            to the given coordinate.
        """
        coord = Vector(coord)
        vs = self.get_verts(coord)
        if len(vs) >= 1:
            # We have hash matches!
            d = 10000000000000.0
            i = -1
            for vi in vs:
                dp = (self.mesh.vertices[vi].co - coord).length
                if dp < d:
                    d = dp
                    i = vi
        else:
            # If all else fails, brute-force.
            d = 10000000000000.0
            i = -1
            for v in self.mesh.vertices:
                dp = (v.co - coord).length
                if dp < d:
                    d = dp
                    i = v.index
        return i


class MeshConnectivityHash:
    """ A hash of mesh connectivity, for easy lookups of mesh elements
        based on relationship.
    """
    def __init__(self, mesh):
        self.mesh = mesh
        self.vhash = {}

        for e in mesh.edges:
            v1 = e.vertices[0]
            v2 = e.vertices[1]

            if v1 not in self.vhash:
                self.vhash[v1] = {"v": [v2], "e": [e.index]}
            else:
                if "v" not in self.vhash[v1]:
                    self.vhash[v1]["v"] = [v2]
                else:
                    self.vhash[v1]["v"] += [v2]
                if "e" not in self.vhash[v1]:
                    self.vhash[v1]["e"] = [e.index]
                else:
                    self.vhash[v1]["e"] += [e.index]

            if v2 not in self.vhash:
                self.vhash[v2] = {"v": [v1], "e": [e.index]}
            else:
                if "v" not in self.vhash[v2]:
                    self.vhash[v2]["v"] = [v1]
                else:
                    self.vhash[v2]["v"] += [v1]
                if "e" not in self.vhash[v1]:
                    self.vhash[v2]["e"] = [e.index]
                else:
                    self.vhash[v2]["e"] += [e.index]

    def v_connected_verts(self, vi):
        """ Returns a list of vertex indices directly connected
            to the given vertex index.
        """
        if vi in self.vhash:
            return self.vhash[vi]["v"]
        else:
            return []


def list_of_affecting_bone_names(obj):
    """ Returns a list of bone names that are deforming
        the mesh through armature modifiers.
    """
    bone_names_d = {}
    for modifier in obj.modifiers:
        if modifier.type == 'ARMATURE':
            arm = modifier.object
            if arm:
                for bone in arm.data.bones:
                    if bone.use_deform:
                        bone_names_d[bone.name] = True
    bone_names = []
    for key in bone_names_d:
        bone_names += [key]
    bone_names.sort()
    return bone_names


def vgroup_mirror_mapping(obj):
    """ Returns a dictionary with key/value pairs corresponding to
        how vgroups should map to each other when mirroring.
        For example, one pair might be {"arm.L": "arm.R"}
    """
    mapping = {}
    for vg in obj.vertex_groups:
        va = vg.name
        if vg.name.endswith(".L"):
            va = vg.name[:-2] + ".R"
        elif vg.name.endswith(".R"):
            va = vg.name[:-2] + ".L"
        elif vg.name.endswith(".r"):
            va = vg.name[:-2] + ".l"
        elif vg.name.endswith(".r"):
            va = vg.name[:-2] + ".l"

        if va in obj.vertex_groups:
            mapping[vg.name] = va
        else:
            mapping[vg.name] = vg.name
    return mapping


def sort_vertices_by_connectivity(vis, mchash):
    """ Sorts the list of vertices based on connectivity.
        If the sorting is abmiguous or if vertices are disconnected,
        it may drop vertices from the list.
        vis: a list of vertex indices
    """
    start_vi = -1
    for vi in vis:
        c = 0
        vcs = mchash.v_connected_verts(vi)
        for vci in vcs:
            if vci in vis:
                c += 1
        if c == 1:
            start_vi = vi
    if start_vi == -1:
        return []

    vlist = [start_vi]
    current_vi = start_vi
    for i in range(len(vis)):
        vcs = mchash.v_connected_verts(current_vi)
        for vci in vcs:
            if vci in vis and vci not in vlist:
                vlist += [vci]
                current_vi = vci
                break

    return vlist


def selected_vertex_indices(obj):
    """ Returns a list of selected vertex indices.
        Assumes you are _not_ in edit mode.
    """
    indices = []
    for v in obj.data.vertices:
        if v.select and not v.hide:
            indices += [v.index]
    return indices


def order_vertex_groups(obj, vi):
    """ Orders the vertex's vertex groups the same as in the master vertex
        group list.
    """
    vgroups = {}
    vertex = obj.data.vertices[vi]
    for group in vertex.groups:
        vgroups[group.group] = group.weight

    for group in vgroups:
        obj.vertex_groups[group].remove([vi])

    for group in obj.vertex_groups:
        if group.index in vgroups:
            group.add([vi], vgroups[group.index], 'REPLACE')


def remove_zerod_vertex_groups(obj, vi):
    """ Removes vertex groups from a vertex that are zero weight.
    """
    vgroups = {}
    vertex = obj.data.vertices[vi]
    for group in vertex.groups:
        vgroups[group.group] = group.weight

    for group in vgroups:
        if abs(vgroups[group]) < 0.0000001:
            obj.vertex_groups[group].remove([vi])


def normalize_bone_weights(obj, vi, vgmap):
    """ Normalizes the bone weights of the given vertex index.
    """
    vgroups = {}
    vertex = obj.data.vertices[vi]

    # Get weights
    for group in vertex.groups:
        gname = obj.vertex_groups[group.group].name
        if gname in vgmap:
            vgroups[gname] = group.weight

    # Normalize weights
    tot = 0.0
    for key in vgroups:
        tot += vgroups[key]
    if tot == 0:
        return
    for key in vgroups:
        vgroups[key] /= tot

    # Set weights to their new values
    for gname in vgroups:
        obj.vertex_groups[gname].add([vi], vgroups[gname], 'REPLACE')


def select_mirror_vertex(obj, vi, vchash):
    """ Selects the x-mirror vertex from the given vertex index.
    """
    coords = obj.data.vertices[vi].co.copy()
    coords[0] = coords[0] * (-1)

    i = vchash.closest_vert(coords)
    if i and vi != i:
        obj.data.vertices[i].select = True


def mirror_selected_vertex_groups(obj, vi, vchash, vgmap):
    """ Copies a vertex's bone vertex groups over to its x-mirror vertex.
        It will mirror vertex group names where appropriate, too.
    """
    v1 = obj.data.vertices[vi]
    vg = obj.vertex_groups

    # Get the mirror vertex index
    mcoords = obj.data.vertices[vi].co.copy()
    mcoords[0] = mcoords[0] * (-1)
    i = vchash.closest_vert(mcoords)

    # Mirror!
    if i:
        if vi != i:
            # Clear the mirror vertex of bone vertex groups
            for vgname in vgmap:
                vg[vgname].remove([i])
            for g in v1.groups:
                gname = vg[g.group].name
                if gname in vgmap:
                    gname = vgmap[gname]
                    vg[gname].add([i], g.weight, 'REPLACE')
        else:
            # Don't touch center vertices
            pass


def mirror_bone_weights(obj, vchash, vgmap, side='LEFT'):
    """ Mirrors all the bone weights of a mesh over the x-axis.
        'LEFT' copies from left to right.
        'RIGHT' copies from right to left
    """
    vg = obj.vertex_groups

    def mirror_side(coord):
        """ Determines if the coordinate is on the proper side
            to receive mirroring.
        """
        if side == 'LEFT':
            return coord[0] < 0
        else:
            return coord[0] > 0

    def mirror_name(gname, receive=False):
        """ Determines if the given vertex group name is a mirror name.
        """
        if gname in vgmap:
            if gname != vgmap[gname]:
                if (side == 'LEFT' and receive) or (side == 'RIGHT' and not receive):
                    return gname.endswith(".r") or gname.endswith(".R")
                else:
                    return gname.endswith(".l") or gname.endswith(".L")
        return False

    for v in obj.data.vertices:
        mco = ((v.co[0] * -1), v.co[1], v.co[2])
        mvi = vchash.closest_vert(mco)
        if mvi == v.index:
            # Center vertex
            for gname in vgmap:
                if mirror_name(gname, receive=True):
                    vg[gname].remove([v.index])
            for g in v.groups:
                gname = vg[g.group].name
                if mirror_name(gname, receive=False):
                    mgname = vgmap[gname]
                    vg[mgname].add([v.index], g.weight, 'REPLACE')
        elif mirror_side(v.co):
            # Vertex to receive mirrored weights
            v_from = obj.data.vertices[mvi]
            v_to = v
            for vgname in vgmap:
                vg[vgname].remove([v_to.index])
            for g in v_from.groups:
                gname = vg[g.group].name
                if gname in vgmap:
                    gname = vgmap[gname]
                    vg[gname].add([v_to.index], g.weight, 'REPLACE')


def get_vertex_groups(obj, vi):
    """ Returns a dict of vertex groups and corresponding weights on the given
        vertex index.
    """
    groups = {}
    v = obj.data.vertices[vi]
    for g in v.groups:
        gname = obj.vertex_groups[g.group].name
        groups[gname] = g.weight
    return groups


def interpolate_bone_weights(obj, vis, vgmap, use_distance=False):
    """ Interpolates bone weights across the list of vertex indices.
    """
    # Create a dict of weight end-point pairs, to interpolate with.
    temp1 = get_vertex_groups(obj, vis[0])
    temp2 = get_vertex_groups(obj, vis[-1])
    for key in temp1:
        if key not in temp2:
            temp2[key] = 0.0
    for key in temp2:
        if key not in temp1:
            temp1[key] = 0.0
    weights = {}
    for key in vgmap:
        if key in temp1:
            weights[key] = (temp1[key], temp2[key])

    # Middle vertices
    mvis = vis[1:][:-1]

    # Remove bone weights from middle vertices
    for gname in vgmap:
        obj.vertex_groups[gname].remove(mvis)

    # Construct alphas for the middle verts
    if not use_distance:
        alphas = [(1.0 / (len(mvis) + 1)) * x for x in range(1, len(mvis) + 1)]
    else:
        dists = []
        for i in range(len(vis) - 1):
            v1 = obj.data.vertices[vis[i]]
            v2 = obj.data.vertices[vis[i + 1]]
            dists += [(v1.co - v2.co).length]
        tot = sum(dists)
        alphas = []
        for i in range(len(dists)):
            alphas += [sum(dists[:i + 1]) / tot]
        alphas = alphas[:-1]

    # Interpolate!
    for vi, a in zip(mvis, alphas):
        for gname in weights:
            w1 = weights[gname][0]
            w2 = weights[gname][1]
            w = w1 + ((w2 - w1) * a)
            obj.vertex_groups[gname].add([vi], w, 'REPLACE')


def select_vertices_with_num_vg(obj, n, and_greater_than=False, and_less_than=False):
    """ Selects the vertices in the mesh with a number of vertex groups
        equal to n.  If and_greater_than is true, it will also select
        vertices with a number of vertex groups greater than n.  And
        respectively less than for and_less_than.
    """
    bm = bmesh.from_edit_mesh(obj.data)

    for vertex in obj.data.vertices:
        if len(vertex.groups) == n:
            bm.verts[vertex.index].select_set(True)
        if and_greater_than and len(vertex.groups) > n:
            bm.verts[vertex.index].select_set(True)
        if and_less_than and len(vertex.groups) < n:
            bm.verts[vertex.index].select_set(True)

    bm.select_flush(True)
    bm.select_flush(False)


def remove_excess_vertex_groups_from_verts(obj, n):
    """ Makes sure all vertices on a mesh have n or fewer vertex groups.
        If the number of vertex groups on a vertex are greater than n, then
        the lowest-weight groups are removed to bring it down to n.
    """
    for vertex in obj.data.vertices:
        if len(vertex.groups) > n:
            # Construct list of weight/index tuples for the vertex groups
            # on the vertex, sorted by weight.
            vg_list = []
            for vg in vertex.groups:
                vg_list += [(vg.weight, vg.group)]
            vg_list.sort()

            # Truncate the list to just the groups to be
            # removed from the vertex.
            diff = len(vertex.groups) - n
            vg_list = vg_list[:diff]

            # Remove the vertex groups
            for (w, g) in vg_list:
                obj.vertex_groups[g].remove([vertex.index])


class MirrorBoneWeights(bpy.types.Operator):
    """ X-axis-mirrors the selected vertices' bone weights over to their
        mirror vertex counterparts.  Based on vertex position, not topology.
        Only bone weights are mirrored.  Other vertex groups are left alone.
    """
    bl_idname = "object.mirror_bone_weights"
    bl_label = "Mirror Bone Weights"

    @classmethod
    def poll(cls, context):
        return context.active_object != None

    def execute(self, context):
        obj = context.active_object
        bpy.ops.object.mode_set(mode='OBJECT')

        bones = list_of_affecting_bone_names(obj)
        vchash = VCoordHash(obj.data)
        vg = vgroup_mirror_mapping(obj)

        vgmap = {}
        for b in bones:
            if b in vg:
                vgmap[b] = vg[b]

        mirror_bone_weights(obj, vchash, vgmap, side='LEFT')

        bpy.ops.object.mode_set(mode='EDIT')
        return {'FINISHED'}


class SortVertexWeights(bpy.types.Operator):
    """ Sorts the selected vertices' vertex weights into the same order as
        the vertex group listing in the mesh properties panel.
    """
    bl_idname = "object.sort_vertex_weights"
    bl_label = "Sort Vertex Weights"

    @classmethod
    def poll(cls, context):
        return context.active_object != None

    def execute(self, context):
        obj = context.active_object
        bpy.ops.object.mode_set(mode='OBJECT')

        for vi in selected_vertex_indices(obj):
            order_vertex_groups(obj, vi)

        bpy.ops.object.mode_set(mode='EDIT')
        return {'FINISHED'}


class InterpolateBoneWeights(bpy.types.Operator):
    """ Interpolates bone weights along a line of selected vertices.
    """
    bl_idname = "object.interpolate_bone_weights"
    bl_label = "Interpolate Bone Weights"

    @classmethod
    def poll(cls, context):
        return context.active_object != None

    def execute(self, context):
        obj = context.active_object
        bones = list_of_affecting_bone_names(obj)
        bpy.ops.object.mode_set(mode='OBJECT')
        mchash = MeshConnectivityHash(obj.data)
        vg = vgroup_mirror_mapping(obj)
        vgmap = {}
        for b in bones:
            if b in vg:
                vgmap[b] = vg[b]

        vis = []
        for vi in selected_vertex_indices(obj):
            vis += [vi]
        l = len(vis)
        if l >= 3:
            vis = sort_vertices_by_connectivity(vis, mchash)
            if len(vis) == l:
                interpolate_bone_weights(obj, vis, vgmap, use_distance=True)

        bpy.ops.object.mode_set(mode='EDIT')
        return {'FINISHED'}


class RemoveZeroWeights(bpy.types.Operator):
    """ Removes vertices from groups for which they have zero weight.
    """
    bl_idname = "object.remove_zero_weights"
    bl_label = "Remove Zero Weights"

    @classmethod
    def poll(cls, context):
        return context.active_object != None

    def execute(self, context):
        obj = context.active_object
        bpy.ops.object.mode_set(mode='OBJECT')

        for vi in selected_vertex_indices(obj):
            remove_zerod_vertex_groups(obj, vi)

        bpy.ops.object.mode_set(mode='EDIT')
        return {'FINISHED'}


class NormalizeBoneWeights(bpy.types.Operator):
    """ Normalizes bone weights in selected vertices.
    """
    bl_idname = "object.normalize_bone_weights"
    bl_label = "Normalize Bone Weights"

    @classmethod
    def poll(cls, context):
        return context.active_object != None

    def execute(self, context):
        obj = context.active_object
        bones = list_of_affecting_bone_names(obj)
        bpy.ops.object.mode_set(mode='OBJECT')
        vg = vgroup_mirror_mapping(obj)
        vgmap = {}
        for b in bones:
            if b in vg:
                vgmap[b] = vg[b]

        for vi in selected_vertex_indices(obj):
            normalize_bone_weights(obj, vi, vgmap)

        bpy.ops.object.mode_set(mode='EDIT')
        return {'FINISHED'}


class SelectMirrorVertices(bpy.types.Operator):
    """ Selects the vertices corresponding the the x-mirrored positions of
        the currently selected vertices.
    """
    bl_idname = "object.select_mirror_vertices"
    bl_label = "Select Mirror Vertices"

    @classmethod
    def poll(cls, context):
        return context.active_object != None

    def execute(self, context):
        obj = context.active_object
        bpy.ops.object.mode_set(mode='OBJECT')
        vchash = VCoordHash(obj.data)

        for vi in selected_vertex_indices(obj):
            select_mirror_vertex(obj, vi, vchash)

        bpy.ops.object.mode_set(mode='EDIT')
        return {'FINISHED'}


class SelectVertsFromNumVG(bpy.types.Operator):
    """ Selects vertices based on the number of vertex groups they have.
    """
    bl_idname = "mesh.select_verts_vg_count"
    bl_options = {'REGISTER', 'UNDO'}
    bl_label = "Select Vertices From Vertex Group Count"

    # Operator Parameters
    num_vg = bpy.props.IntProperty(
               name="num_vg",
               default=5, min=0, soft_min=0, soft_max=20,
               description="Number of vertex groups",
             )
    mode_items = [('EQUAL_OR_LESS',    "Equal or Less Than",    "", 0),
                  ('EQUAL',            "Equal",                 "", 1),
                  ('EQUAL_OR_GREATER', "Equal or Greater Than", "", 2)]
    mode = bpy.props.EnumProperty(
             items=mode_items,
             name="Mode",
             default='EQUAL_OR_GREATER'
           )

    # Operator Methods
    @classmethod
    def poll(cls, context):
        return context.active_object != None and context.mode == 'EDIT_MESH'

    def execute(self, context):
        obj = context.active_object

        bpy.ops.object.mode_set(mode='OBJECT')
        bpy.ops.object.mode_set(mode='EDIT')
        if self.mode == 'EQUAL':
            select_vertices_with_num_vg(obj, self.num_vg, False, False)
        elif self.mode == 'EQUAL_OR_LESS':
            select_vertices_with_num_vg(obj, self.num_vg, False, True)
        elif self.mode == 'EQUAL_OR_GREATER':
            select_vertices_with_num_vg(obj, self.num_vg, True, False)

        return {'FINISHED'}


class RemoveExcessVGFromVerts(bpy.types.Operator):
    """ Reduces the number of vertex groups on each vertex down to a
        user-specified maximum.  Only the lowest-weight groups are removed
        to bring them down to the maximum.
    """
    bl_idname = "mesh.remove_excess_vg_from_verts"
    bl_options = {'REGISTER', 'UNDO'}
    bl_label = "Remove Excess Vertex Groups From Verts"

    # Operator Parameters
    max_vg = bpy.props.IntProperty(
               name="max_vg",
               default=4, min=0, soft_min=0, soft_max=20,
               description="Max number of vertex groups to allow",
             )

    # Operator Methods
    @classmethod
    def poll(cls, context):
        return context.active_object != None and context.mode == 'EDIT_MESH'

    def execute(self, context):
        obj = context.active_object

        bpy.ops.object.mode_set(mode='OBJECT')
        remove_excess_vertex_groups_from_verts(obj, self.max_vg)
        bpy.ops.object.mode_set(mode='EDIT')

        return {'FINISHED'}


class VertexGroupStatistics(bpy.types.Panel):
    bl_space_type = 'VIEW_3D'
    bl_region_type = 'UI'
    bl_label = "Vertex Group Statistics"
    bl_idname = "PT_vertex_group_statistics"

    @classmethod
    def poll(self, context):
        return context.active_object != None and context.mode == 'EDIT_MESH'

    def draw(self, context):
        # Calculate various statistics
        bm = bmesh.from_edit_mesh(context.active_object.data)

        selected_vert_count = 0

        max_weight = 0.0
        min_weight = 1.0
        avg_weight = 0.0

        max_groups = 0
        min_groups = 99999
        avg_groups = 0

        for bv in bm.verts:
            if bv.select:
                selected_vert_count += 1
                v = context.active_object.data.vertices[bv.index]
                ng = len(v.groups)

                max_groups = max(max_groups, ng)
                min_groups = min(min_groups, ng)
                avg_groups += ng

        avg_groups /= selected_vert_count

        layout = self.layout

        row = layout.row()
        row.label("In selected vertices:")

        row = layout.row()
        row.label("Max groups: " + str(max_groups))
        row = layout.row()
        row.label("Min groups: " + str(min_groups))
        row = layout.row()
        row.label("Average groups: " + str(avg_groups))


def register():
    bpy.utils.register_class(SelectVertsFromNumVG)
    bpy.utils.register_class(RemoveExcessVGFromVerts)
    bpy.utils.register_class(VertexGroupStatistics)
    bpy.utils.register_class(MirrorBoneWeights)
    bpy.utils.register_class(SortVertexWeights)
    bpy.utils.register_class(InterpolateBoneWeights)
    bpy.utils.register_class(NormalizeBoneWeights)
    bpy.utils.register_class(SelectMirrorVertices)
    bpy.utils.register_class(RemoveZeroWeights)


def unregister():
    bpy.utils.unregister_class(SelectVertsFromNumVG)
    bpy.utils.unregister_class(RemoveExcessVGFromVerts)
    bpy.utils.unregister_class(VertexGroupStatistics)
    bpy.utils.unregister_class(MirrorBoneWeights)
    bpy.utils.unregister_class(SortVertexWeights)
    bpy.utils.unregister_class(InterpolateBoneWeights)
    bpy.utils.unregister_class(NormalizeBoneWeights)
    bpy.utils.unregister_class(SelectMirrorVertices)
    bpy.utils.unregister_class(RemoveZeroWeights)


if __name__ == "__main__":
    register()
