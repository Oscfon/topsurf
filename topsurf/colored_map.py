r"""
Vertex colored combinatorial maps

The main class in this module is :class:`ColoredOrientedMap`.
"""
# ****************************************************************************
#  This file is part of topsurf
#
#       Copyright (C) 2026 Juliette Schabanel
#
#  This program is free software; you can redistribute it and/or
#  modify it under the terms of the GNU General Public License
#  as published by the Free Software Foundation; either version 2
#  of the License, or (at your option) any later version.
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  You should have received a copy of the GNU General Public License
#  along with this program; if not, write to the Free Software
#  Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.
# ****************************************************************************

from topsurf import OrientedMap
from topsurf.permutation import (perm_cycles, perm_cycles_to_string, perm_orbit)

class ColoredOrientedMap(OrientedMap):
    r"""
    INPUT:
        - ``ecolors`` -- ``None`` or a list or dictionnary of data so that ecolors[i] is the data associated to the edge (2i,2i+1).
        - ``vcolors`` -- ``None`` or a dictionnary of data so that vcolors[h] is the data associated to the vertex incident to half edge half, and only one half edge incident to a given vertex is present.
    """
    def __init__(self, vp=None, fp=None, vcolors=None, ecolors=None, mutable=False, check=True):
        OrientedMap.__init__(self, vp, fp, mutable, check)

        if ecolors is None:
            self._edge_colors = [None]*(len(self._vp)//2)
        else:
            if isinstance(ecolors, list):
                if check:
                    if len(ecolors) != len(self._vp)//2:
                        raise ValueError("The number of colors doesn't match the number of edges")
                self._edge_colors = ecolors
            elif isinstance(ecolors, dict):
                self._edge_colors = [ecolors.get(e) for e in range(len(self._vp)//2)]
            else:
                raise TypeError("ecolors should be a list or dictionnary")

        if len(self._vp) == 0:
            if vcolors is None:
                self._vertex_colors = {-1: None}
            else:
                self._vertex_colors = {-1: vcolors.get(-1)}
        else:
            self._vertex_colors = {}
            for v in self.vertices():
                col = None
                for h in v:
                    hcol = vcolors.get(h)
                    if hcol is not None:
                        if col is None:
                            col = hcol
                        else:
                            raise ValueError(f"multiple values for {v} in vcolors")
                self._vertex_colors[v[0]] = col

    def __str__ (self):
        vp_cycles = perm_cycles(self._vp)
        vertices = perm_cycles_to_string(vp_cycles, edge_like=True)
        fp_cycles = perm_cycles(self._fp)
        faces = perm_cycles_to_string(fp_cycles, edge_like=True)
        return f"ColoredOrientedMap(\"{vertices}\", \"{faces}\", \n edge colors: {self._edge_colors}, \n vertex colors: {self._vertex_colors})"

    __repr__=__str__

    def copy(self, mutable=None):
        r"""
        Return a copy of this oriented map.

        If ``mutable`` is set to ``True`` or ``False`` return respectively an
        immutable or mutable version of the map.
        """

        if mutable is None:
            mutable = self._mutable

        if not self._mutable and not mutable:
            # avoid copies of immutable objects
            return self

        m = ColoredOrientedMap.__new__(ColoredOrientedMap)
        m._fp = self._fp[:]
        m._vp = self._vp[:]
        m._mutable = mutable
        m._edge_colors = self._edge_colors
        m._vertex_colors = self._vertex_colors
        return m

    def edge_color(self, e):
        r"""
        Return the color of the edge e.
        """
        return self._edge_colors[e]

    def _v_id(self, h):
        r"""
        Return the id of the vertex v incident to h (= half edge of minimal label incident to v).
        """
        h = self._check_half_edge(h)
        return min(perm_orbit(self._vp, h))

    def vertex_color(self, h):
        r"""
        Return the color of the vertex incident to h.
        """
        return self._vertex_colors[self._v_id(h)]

    def set_edge_color(self, h, col, check=True)
        r"""
        Change the color of the half edge h to col
        """
        if check:
            self._assert_mutable()
        self._edge_colors[h//2] = col

    def set_vertex_color(self, h, col, check=True)
        r"""
        Change the color of the vertex incident to half edge h to col
        """
        if check:
            self._assert_mutable()
        self._vertex_colors[self._v_id(h)] = col

    
    def submap(self, edges, mutable=False, check=True):
        r"""
        Return the submap of this constellation induced on ``edges``.

        EXAMPLES::
        """

        unlab_submap = OrientedMap.submap(self, edges, mutable=mutable, check=check)
        relab_submap = OrientedMap.submap(self, edges, relabel=True, mutable=mutable, check=check)
        
        sub_edge_colors = [self.edge_color(e) for e in edges]
        sub_vert_colors = {(v[0] - 2*sum([e not in edges for e in range(v[0]//2)])): self.vertex_color(v[0]) for v in unlab_submap.vertices()}
        
        return ColoredOrientedMap(vp=relab_submap._vp, ecolors=sub_edge_colors, vcolors=sub_vert_colors, mutable=mutable, check=check)
            
    
            






        