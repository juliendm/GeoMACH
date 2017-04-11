"""
B-spline surface-modeling engine (BSE) vec class
"""

# pylint: disable=E1101

from __future__ import division
import numpy

import operator

class BSEvec(object):

    def __init__(self, name, size, ndim, hidden):
        self.array = numpy.zeros((size, ndim))
        self.name = name
        self.size = size
        self.ndim = ndim
        self._hidden = hidden
        self._file = None
        self._default_var_names = ['v' + str(idim) 
                                   for idim in xrange(self.ndim)]

    def _open_file(self, filename):
        self._file = open(filename, 'w')

    def _write_tec_header(self, title, variables):
        self._file.write('title = ' + title + '\n')
        self._file.write('variables = ')
        for ivar in range(len(variables)):
            self._file.write(variables[ivar] + ',')
        self._file.write('\n')

    def _write(self, text):
        self._file.write(text)

    def _write_line(self, data, label=''):
        self._file.write(label)
        for ind in range(data.shape[0]):
            if data[ind] == data[ind]:
                self._file.write(str(data[ind]) + ' ')
            else:
                self._file.write(str(0.0) + ' ')
        self._file.write('\n')

    def _close_file(self):
        self._file.close()

    def export_tec_scatter(self, filename=None, var_names=None):
        if filename is None:
            filename = self.name + '_scatter.dat'
        if var_names is None:
            var_names = self._default_var_names

        self._open_file(filename)
        self._write_tec_header('BSE output', var_names)
        for ind in xrange(self.size):
            self._write_line(self.array[ind, :])
        self._close_file()



class BSEvecUns(BSEvec):

    pass



class BSEvecStr(BSEvec):

    def __init__(self, name, size, ndim, surf_sizes, hidden,
                 bse=None):
        super(BSEvecStr, self).__init__(name, size, ndim, hidden)

        self._bse = bse
        self.surfs = []
        if surf_sizes is not None:
            ind1, ind2 = 0, 0
            for isurf in xrange(surf_sizes.shape[0]):
                num_u, num_v = surf_sizes[isurf, :]
                ind2 += num_u * num_v
                surf = self.array[ind1:ind2]
                surf = surf.reshape((num_u, num_v, ndim), 
                                    order='F')
                ind1 += num_u * num_v
                self.surfs.append(surf)

    def __call__(self, isurf):
        return self.surfs[isurf]

    def export_tec_str(self, filename=None, var_names=None):
        if filename is None:
            filename = self.name + '_surf.dat'
        if var_names is None:
            var_names = self._default_var_names

        self._open_file(filename)
        self._write_tec_header('BSE output', var_names)
        for isurf in xrange(len(self.surfs)):
            if not self._hidden[isurf]:
                surf = self.surfs[isurf]
                num_u, num_v = surf.shape[:2]
                self._write('zone i='+str(num_u) + \
                            ', j=' + str(num_v) + \
                            ', DATAPACKING=POINT\n')
                for ind_v in xrange(num_v):
                    for ind_u in xrange(num_u):
                        self._write_line(surf[ind_u, ind_v, :])
        self._close_file()

    def export_STL(self, filename=None, var_names=None):
        if filename is None:
            filename = self.name + '.stl'
        if var_names is None:
            var_names = self._default_var_names

        self._open_file(filename)
        self._write('solid model\n')

#        for surf in self.surfs:
        for isurf in xrange(len(self.surfs)):
            if not self._hidden[isurf]:
                surf = self.surfs[isurf]

                num_u, num_v = surf.shape[:2]
                for ind_v in xrange(num_v - 1):
                    for ind_u in xrange(num_u - 1):
                        pt1 = surf[ind_u, ind_v+1, :]
                        pt2 = surf[ind_u, ind_v, :]
                        pt3 = surf[ind_u+1, ind_v, :]
                        
                        for pts in [[surf[ind_u, ind_v+1, :],
                                     surf[ind_u, ind_v, :],
                                     surf[ind_u+1, ind_v, :]],
                                    [surf[ind_u+1, ind_v, :],
                                     surf[ind_u+1, ind_v+1, :],
                                     surf[ind_u, ind_v+1, :]]]:
                            nor = numpy.cross(pts[2] - pts[1], 
                                              pts[0] - pts[1])
                            nor /= numpy.linalg.norm(nor)
                            self._write_line(nor, 'facet normal ')
                            self._write('outer loop\n')
                            self._write_line(pts[0], 'vertex ')
                            self._write_line(pts[1], 'vertex ')
                            self._write_line(pts[2], 'vertex ')
                            self._write('endloop\n')
                            self._write('endfacet\n')
        self._write('endsolid model')
        self._close_file()

    def export_MESH(self, filename=None):
        if filename is None:
            filename = 'surface_fluid.mesh'

        fluid = open(filename, 'w')

        nElem = 0

        pts = []
        pts_filtered = []
        corresp = []
        lumped = []

        index = 0
        for isurf in xrange(len(self.surfs)):
            if not self._hidden[isurf]:
                surf = self.surfs[isurf]
                num_u, num_v = surf.shape[:2]
                nElem += (num_u-1)*(num_v-1)
                for ind_v in xrange(num_v):
                    for ind_u in xrange(num_u):
                        pt = surf[ind_u, ind_v, :]
                        pts.append([round(pt[0],9),round(pt[1],9),round(pt[2],9),index])
                        index += 1
        pts_sorted = sorted(pts, key = operator.itemgetter(0, 1, 2))

        eps = 1e-8
        index = pts_sorted[0][3]
        corresp.append([pts_sorted[0][3], index])
        for ipt in range(1,len(pts_sorted)):
            if (abs(pts_sorted[ipt][0]-pts_sorted[ipt-1][0]) > eps or abs(pts_sorted[ipt][1]-pts_sorted[ipt-1][1]) > eps or abs(pts_sorted[ipt][2]-pts_sorted[ipt-1][2]) > eps):
                index = pts_sorted[ipt][3]
            corresp.append([pts_sorted[ipt][3], index])
        corresp.sort(key = operator.itemgetter(0))

        pts_filtered.append(pts[0])
        index = 0
        for ipt in range(1,len(pts)):
            lumped.append(index)
            if (corresp[ipt][0] == corresp[ipt][1]):
                pts_filtered.append(pts[ipt])
                index += 1

        nVertex = len(pts_filtered)

        string = '\nMeshVersionFormatted\n2\n\nDimension\n3\n\nVertices\n' + str(nVertex) + '\n\n'
        fluid.write(string)

        for pt in pts_filtered:
            y = pt[2]
            if (abs(y) < 0.001): y = 0.0
            string = str(pt[0]) + " " + str(y) + " " + str(pt[1]) + " 0\n" # ROTATE FOR CFD
            fluid.write(string)

        remove = [] # [455] #[452, 364, 367, 370, 373] #[450, 96, 456, 454, 452, 373, 372, 371, 370, 369, 368, 367, 366, 365, 364, 363, 362] #[96, 92, 451, 450, 453, 455, 457]

        # wing_bottom = [388,392,396,389,393,397,111,110,109,113,112,440,115,114,442,444,116,117,119,118,446,448,120,121,123,122,450,452,124,125,127,126,454,456,128,129,131,130,458,114,242,243,244,253,252,251,260,261,262,271,270,269,278,279,280,289,288,287,298,297,296,305,306,307,316,315,314,323,324,325,334,333,332,341,342,343,423,421,419,245,254,263,272,281,290,299,308,317,326,335,344,425,427,345,336,327,318,309,300,291,282,273,264,255,246,247,256,265,274,283,292,301,310,319,328,337,346,429,431,347,338,329,320,311,302,293,284,275,266,257,248,249,258,267,276,285,294,303,312,321,330,339,348,433,435,349,340,331,322,313,304,295,286,277,268,259,250]
        # wing_top = [434,432,430,428,426,424,422,420,418,241,240,239,238,237,236,235,234,233,224,225,226,227,228,229,230,231,232,223,222,221,220,219,218,217,216,215,206,207,208,209,210,211,212,213,214,205,204,203,202,201,200,199,198,197,188,189,190,191,192,193,194,195,196,187,186,185,184,183,182,181,180,179,170,171,172,173,174,175,176,177,178,169,168,167,166,165,164,163,162,161,152,153,154,155,156,157,158,159,160,151,150,149,148,147,146,145,144,143,134,135,136,137,138,139,140,141,142]

        fluid.write('\nTriangles\n' + str(2*(nElem - len(remove)*((num_u-1)*(num_v-1))) + 23) + '\n\n')

        base_index = 0
        ref = 0
        for isurf in xrange(len(self.surfs)):
            if not self._hidden[isurf]:
                ref = ref+1;
                surf = self.surfs[isurf]
                num_u, num_v = surf.shape[:2]

                if ref not in remove:
                    for ind_v in xrange(num_v - 1):
                        for ind_u in xrange(num_u - 1):

                            ids = numpy.array([ind_u+(ind_v+1)*num_u, ind_u+ind_v*num_u, (ind_u+1)+ind_v*num_u, (ind_u+1)+(ind_v+1)*num_u]) + base_index

                            coord_0 = pts_filtered[lumped[corresp[ids[0]][1]]]
                            coord_1 = pts_filtered[lumped[corresp[ids[1]][1]]]
                            coord_3 = pts_filtered[lumped[corresp[ids[3]][1]]]

                            # angle 01 - 03
                            D01 = ( (coord_0[0]-coord_1[0])**2.0 + (coord_0[1]-coord_1[1])**2.0 + (coord_0[2]-coord_1[2])**2.0 )**0.5
                            D03 = ( (coord_0[0]-coord_3[0])**2.0 + (coord_0[1]-coord_3[1])**2.0 + (coord_0[2]-coord_3[2])**2.0 )**0.5
                            D13 = ( (coord_1[0]-coord_3[0])**2.0 + (coord_1[1]-coord_3[1])**2.0 + (coord_1[2]-coord_3[2])**2.0 )**0.5

                            # print D01, D03, D13

                            angle = numpy.arccos((D01**2.0 + D03**2.0 - D13**2.0) / (2.0 * D01 * D03))

                            # print angle

                            if (angle >= 0.5*numpy.pi):

                                fluid.write(str(lumped[corresp[ids[0]][1]]+1) + " " + str(lumped[corresp[ids[1]][1]]+1)  + " " + str(lumped[corresp[ids[2]][1]]+1) + " " + str(ref) + "\n")
                                fluid.write(str(lumped[corresp[ids[2]][1]]+1) + " " + str(lumped[corresp[ids[3]][1]]+1)  + " " + str(lumped[corresp[ids[0]][1]]+1) + " " + str(ref) + "\n")
                            
                            else:

                                fluid.write(str(lumped[corresp[ids[1]][1]]+1) + " " + str(lumped[corresp[ids[2]][1]]+1)  + " " + str(lumped[corresp[ids[3]][1]]+1) + " " + str(ref) + "\n")
                                fluid.write(str(lumped[corresp[ids[3]][1]]+1) + " " + str(lumped[corresp[ids[0]][1]]+1)  + " " + str(lumped[corresp[ids[1]][1]]+1) + " " + str(ref) + "\n")

                base_index += num_u*num_v

        #fluid.write('30837 30836 31674 1000\n30836 31675 31674 1000\n30836 30835 31675 1000\n30835 31676 31675 1000\n30835 30834 31676 1000\n30834 31677 31676 1000\n30834 30833 31677 1000\n30833 31678 31677 1000\n30833 30832 31678 1000\n30832 31679 31678 1000\n30832 30831 31679 1000\n30831 31680 31679 1000\n30831 30830 31680 1000\n30830 31681 31680 1000\n30830 30829 31681 1000\n30829 31682 31681 1000\n30829 30756 31682 1000\n30756 31683 31682 1000\n30756 30755 31683 1000\n30755 31756 31683 1000\n30755 30754 31756 1000\n30754 31757 31756 1000\n30754 30753 31757 1000\n30753 31758 31757 1000\n30753 30752 31758 1000\n30752 31759 31758 1000\n30752 30751 31759 1000\n30751 31760 31759 1000\n30751 30750 31760 1000\n30750 31761 31760 1000\n30750 30749 31761 1000\n30749 31762 31761 1000\n30749 30748 31762 1000\n30748 31763 31762 1000\n30748 30675 31763 1000\n30675 31764 31763 1000\n30675 30674 31764 1000\n30674 31829 31764 1000\n30674 30673 31829 1000\n30673 31830 31829 1000\n30673 30672 31830 1000\n30672 31831 31830 1000\n30672 30671 31831 1000\n30671 31832 31831 1000\n30671 30670 31832 1000\n30670 31833 31832 1000\n30670 30669 31833 1000\n30669 31834 31833 1000\n30669 30668 31834 1000\n30668 31835 31834 1000\n30668 30667 31835 1000\n30667 31836 31835 1000\n30667 30666 31836 1000\n')
        
        fluid.write('9551 9337 9338 1000\n9551 9552 9337 1000\n9552 9336 9337 1000\n9552 9553 9336 1000\n9553 9335 9336 1000\n9553 9554 9335 1000\n9554 9318 9335 1000\n9554 9555 9318 1000\n9555 9317 9318 1000\n9555 9572 9317 1000\n9572 9316 9317 1000\n9572 9573 9316 1000\n9573 9315 9316 1000\n9573 9574 9315 1000\n9574 9298 9315 1000\n9574 9575 9298 1000\n9575 9297 9298 1000\n9575 9588 9297 1000\n9588 9296 9297 1000\n9588 9589 9296 1000\n9589 9295 9296 1000\n9589 9590 9295 1000\n9590 9294 9295 1000')

        fluid.write('\nEnd\n')

        fluid.close()

    def export_IGES(self, filename=None):
        ks = []
        ms = []
        ds = [[],[]]
        Cs = []

        for isurf in xrange(len(self.surfs)):
            if not self._hidden[isurf]:
                ku = self._bse.get_bspline_option('order', isurf, 'u')
                kv = self._bse.get_bspline_option('order', isurf, 'v')
                mu = self._bse.get_bspline_option('num_cp', isurf, 'u')
                mv = self._bse.get_bspline_option('num_cp', isurf, 'v')
                du = numpy.zeros(ku + mu)
                dv = numpy.zeros(kv + mv)
                du[mu:] = 1.0
                dv[mv:] = 1.0
                du[ku-1:mu+1] = numpy.linspace(0, 1, mu-ku+2)
                dv[kv-1:mv+1] = numpy.linspace(0, 1, mv-kv+2)

                ks.append([ku,kv])
                ms.append([mu,mv])
                ds[0].append(du)
                ds[1].append(dv)
                Cs.append(self.surfs[isurf])
        ks = numpy.array(ks)
        ms = numpy.array(ms)

        def write(val, dirID, parID, field, last=False):
            if last:
                self._write('%20.12e;' %(val.real))
            else:
                self._write('%20.12e,' %(val.real))
            field += 1
            if field==3:
                field = 0
                self._write('%9i' %(dirID))
                self._write('P')
                self._write('%7i\n' %(parID))
                parID += 1
            return parID, field

        if filename is None:
            filename = self.name + '.igs'

        self._open_file(filename)
        self._write('                                                                        S      1\n')
        self._write('1H,,1H;,4HSLOT,37H$1$DUA2:[IGESLIB.BDRAFT.B2I]SLOT.IGS;,                G      1\n')
        self._write('17HBravo3 BravoDRAFT,31HBravo3->IGES V3.002 (02-Oct-87),32,38,6,38,15,  G      2\n')
        self._write('4HSLOT,1.,1,4HINCH,8,0.08,13H871006.192927,1.E-06,6.,                   G      3\n')
        self._write('31HD. A. Harrod, Tel. 313/995-6333,24HAPPLICON - Ann Arbor, MI,4,0;     G      4\n')

        dirID = 1
        parID = 1
        for s in range(ks.shape[0]):
            numFields = 4 + ds[0][s].shape[0] + ds[1][s].shape[0] + 4*ms[s,0]*ms[s,1]
            numLines = 2 + numpy.ceil(numFields/3.0)
            for val in [128, parID, 0, 0, 1, 0, 0, 0]:
                self._write('%8i' %(val))
            self._write('00000001')
            self._write('D')
            self._write('%7i\n' %(dirID))
            dirID += 1
            for val in [128, 0, 2, numLines, 0]:
                self._write('%8i' %(val))
            self._write('%32i' %(0))
            self._write('D')
            self._write('%7i\n' %(dirID))
            dirID += 1
            parID += numLines
        nDir = dirID - 1

        dirID = 1
        parID = 1
        for s in range(ks.shape[0]):
            ku = ks[s,0]
            kv = ks[s,1]
            mu = ms[s,0]
            mv = ms[s,1]
            du = ds[0][s]
            dv = ds[1][s]

            for val in [128, mu-1, mv-1, ku-1, kv-1]:
                self._write('%10i,' %(val))  
            self._write('          ')
            self._write('%7i' %(dirID))   
            self._write('P')
            self._write('%7i\n' %(parID))
            parID += 1

            for val in [0, 0, 1, 0, 0]:
                self._write('%10i,' %(val))
            self._write('          ')
            self._write('%7i' %(dirID))   
            self._write('P')
            self._write('%7i\n' %(parID))
            parID += 1

            field = 0
            for i in range(du.shape[0]):
                parID,field = write(du[i], dirID, parID, field)
            for i in range(dv.shape[0]):
                parID,field = write(dv[i], dirID, parID, field)
            for i in range(mu*mv):
                parID,field = write(1.0, dirID, parID, field)
            for j in range(mv):
                for i in range(mu):
                    for k in range(3):
                        parID,field = write(Cs[s][i,j,k], dirID, parID, field)
            parID,field = write(0, dirID, parID, field)
            parID,field = write(1, dirID, parID, field)
            parID,field = write(0, dirID, parID, field)
            parID,field = write(1, dirID, parID, field, last=True)
            if not field==0:
                for i in range(3-field):
                    self._write('%21s' %(' '))
                self._write('%9i' %(dirID))
                self._write('P')
                self._write('%7i\n' %(parID))
                parID += 1

            dirID += 2

        nPar = parID - 1

        self._write('S%7i' %(1))   
        self._write('G%7i' %(4))   
        self._write('D%7i' %(nDir))   
        self._write('P%7i' %(nPar))   
        self._write('%40s' %(' '))   
        self._write('T')
        self._write('%7i\n' %(1))
        self._close_file()
