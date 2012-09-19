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
    "version": (0,1),
    "blender": (2, 6, 3),
    "api": 50557,
    "location": "",
    "wiki_url": "",
    "tracker_url": "",
    "category": "Rigging"}

import bpy
import bmesh




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
        
    
class SelectVertsFromNumVG(bpy.types.Operator):
    """ Selects vertices based on the number of vertex groups they have.
    """
    bl_idname = "mesh.select_verts_vg_count"
    bl_options = {'REGISTER', 'UNDO'}
    bl_label = "Select Vertices From Vertex Group Count"

    # Operator Parameters
    num_vg = bpy.props.IntProperty(
               name="num_vg",
               default=5, min = 0, soft_min = 0, soft_max = 20,
               description="Number of vertex groups",
             )
    mode_items = [('EQUAL_OR_LESS',    "Equal or Less Than",    "", 0),
                  ('EQUAL',            "Equal",                 "", 1),
                  ('EQUAL_OR_GREATER', "Equal or Greater Than", "", 2)]
    mode = bpy.props.EnumProperty(
             items = mode_items,
             name = "Mode",
             default = 'EQUAL_OR_GREATER'
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
            for (w,g) in vg_list:
	            obj.vertex_groups[g].remove([vertex.index])
        
    
class RemoveExcessVGFromVerts(bpy.types.Operator):
    """ Makes sure all vertices on a mesh have n or fewer vertex groups.
        If the number of vertex groups on a vertex are greater than n, then
        the lowest-weight groups are removed to bring it down to n.
    """
    bl_idname = "mesh.remove_excess_vg_from_verts"
    bl_options = {'REGISTER', 'UNDO'}
    bl_label = "Remove Excess Vertex Groups From Verts"

    # Operator Parameters
    max_vg = bpy.props.IntProperty(
               name="max_vg",
               default=4, min = 0, soft_min = 0, soft_max = 20,
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
    
def unregister():
    bpy.utils.unregister_class(SelectVertsFromNumVG)
    bpy.utils.unregister_class(RemoveExcessVGFromVerts)
    bpy.utils.unregister_class(VertexGroupStatistics)


if __name__ == "__main__":
    register()
