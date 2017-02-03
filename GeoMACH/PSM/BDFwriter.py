from __future__ import division
import numpy


def writeBDF(filename, nodes, quads, symm, quad_groups, group_names, group_names_oml,
             new_mem, new_nodes, new_ucoord, new_vcoord):
    f = open(filename, 'w')

    def writeLine(line):
        write(line,r=80)
        write('\n')

    def write(line,l=0,r=0):
        if l is not 0:
            n = l - len(line)
            for i in range(n):
                line = ' ' + line
        if r is not 0:
            n = r - len(line)
            for i in range(n):
                line = line + ' '
        f.write(line)

    unique = numpy.unique(quad_groups)

    new = numpy.zeros(len(quad_groups), int)
    k = 0
    for i in xrange(1+numpy.max(quad_groups)):
        if numpy.prod(i!=quad_groups) == 0:
            new += k * (i==quad_groups)
            k += 1
    quad_groups[:] = new

    writeLine('SOL 200')
    writeLine('TIME 600')
    writeLine('CEND')
    writeLine('SPC = 1')
    writeLine('ANALYSIS = STATICS')
    writeLine('DESOBJ = 1')
    writeLine('DESSUB = 1')
    writeLine('DISPLACEMENT = ALL')
    writeLine('FORCE = ALL')
    writeLine('STRESS = ALL')
    writeLine('SUBCASE 1')
    writeLine('    LOAD = 1')
    writeLine('BEGIN BULK')

    for i in xrange(len(unique)):
        write('$CDSCRPT')
        write(str(i+1), l=16)
        name = group_names[unique[i]]
        write(name+'/'+name[:name.index(':')], l=32)
        write('\n')

    writeLine('$')
    writeLine('$       grid data              0')

    nElem = 0
    used_vertice = numpy.zeros(nodes.shape[0], bool)
    used_elem = numpy.zeros(quads.shape[0], bool)
    for i in range(quads.shape[0]):
        nElem += 1
        used_elem[i] = True
        for k in xrange(4):
            used_vertice[quads[:,k]-1] = True
    names_id = []
    for i in xrange(len(unique)):
        names_id.append(group_names[unique[i]])
    nElem_filtered = 0
    used_vertice_filtered = numpy.zeros(nodes.shape[0], bool)
    used_elem_filtered = numpy.zeros(quads.shape[0], bool)
    for i in range(quads.shape[0]):
        if names_id[quad_groups[i]] in group_names_oml:
            nElem_filtered += 1
            used_elem_filtered[i] = True
            for k in xrange(4):
                used_vertice_filtered[quads[i,k]-1] = True

    nVertices = 0
    node_indices = numpy.zeros(nodes.shape[0], int)
    nVertices_filtered = 0
    node_indices_filtered = numpy.zeros(nodes.shape[0], int)
    for k in xrange(nodes.shape[0]):
        if used_vertice[k]:
            nVertices += 1
            node_indices[k] = nVertices
        if used_vertice_filtered[k]:
            nVertices_filtered += 1
            node_indices_filtered[k] = nVertices_filtered

    for k in range(nodes.shape[0]):
        if used_vertice[k]:
            write('GRID*   ')
            write(str(node_indices[k]),l=16)
            write('0',l=16)
            write('%.8E' % nodes[k,0],l=16)
            write('%.8E' % nodes[k,1],l=16)
            write('*')
            write(str(node_indices[k]),l=7)
            write('\n')
            write('*')
            write(str(node_indices[k]),l=7)
            write('%.8E' % nodes[k,2],l=16)
            write('0',l=16)
            write(' ',l=16)
            write('0',l=16)
            write(' ',l=8)
            write('\n')

    for i in range(quads.shape[0]):
        write('CQUAD4  ')
        write(str(i+1),l=8)
        write(str(quad_groups[i]+1),l=8)
        #write(str(i+1),l=8)
        write(str(node_indices[quads[i,0]-1]),l=8)
        write(str(node_indices[quads[i,1]-1]),l=8)
        write(str(node_indices[quads[i,2]-1]),l=8)
        write(str(node_indices[quads[i,3]-1]),l=8)
        write('\n')

        q = quads[i,:]
        if q[0]==q[1] or q[0]==q[2] or q[0]==q[3] or \
           q[1]==q[2] or q[1]==q[3] or q[2]==q[3]:
            print 'invalid quad', q, group_names[quad_groups[i]], quad_groups[i]

        imem = quad_groups[i]
        uv_selector = imem == new_mem
        for k in range(0):#4):
            inode = quads[i,k] - 1
            node = nodes[inode, :]
            candidates = numpy.array(new_nodes)
            for j in range(3):
                candidates[:,j] -= node[j]
            candidates = candidates ** 2
            candidates = numpy.sum(candidates, axis=1)
            candidates += ~uv_selector * 1e16
            imin = numpy.argmin(candidates)
            write('$PARAM'+str(k)+' ')
            write(' ',l=16)
            write('%.8E' % new_ucoord[imin],l=16)
            write('%.8E' % new_vcoord[imin],l=16)
            #write(' ',l=24)
            write('\n')

    min_x_index = 0
    min_x = 1.0e20
    max_x_index = 0
    max_x = -1.0e20
    for i in range(nodes.shape[0]):
        if symm[i] and used_vertice[i]:
            if nodes[i,0] > max_x:
                max_x = nodes[i,0]
                max_x_index = node_indices[i]
            if nodes[i,0] < min_x:
                min_x = nodes[i,0]
                min_x_index = node_indices[i]

    for i in range(nodes.shape[0]):
        if symm[i] and used_vertice[i]: # and used_vertice_filtered[i]: ##################################
            write('SPC     ')
            write('1',l=8)
            write(str(node_indices[i]),l=8)
            if used_vertice_filtered[i]: #node_indices[i] in [min_x_index,max_x_index]:
                write('  123456')
            else:
                write('       3')
            write('     0.0')
            #write(' ',l=40)
            write('\n')

    writeLine('END BULK')

    f.close()

    writeMESH(filename.split('.bdf')[0] + '_surface.mesh', nodes, quads, quad_groups, node_indices_filtered, nVertices_filtered, nElem_filtered, used_vertice_filtered, used_elem_filtered, used_elem_filtered)
    writeMESH(filename.split('.bdf')[0] + '.mesh', nodes, quads, quad_groups, node_indices, nVertices, nElem, used_vertice, used_elem, used_elem_filtered)

def writeMESH(filename, nodes, quads, quad_groups, node_indices, nVertices, nElem, used_vertice, used_elem, used_elem_filtered):

    struct = open(filename, 'w')

    struct.write('\nMeshVersionFormatted\n2\n\nDimension\n3\n\nVertices\n' + str(nVertices) + '\n\n')

    for k in range(nodes.shape[0]):
        if used_vertice[k]:
            struct.write(str(nodes[k,0]) + " " + str(nodes[k,2]) + " " + str(nodes[k,1]) + " "+ str(node_indices[k]) +"\n")

    struct.write('\nQuadrilaterals\n' + str(nElem) + '\n\n')

    for k in range(quads.shape[0]):
        if used_elem[k]:
            if used_elem_filtered[k]:
                struct.write(str(node_indices[quads[k,0]-1]) + " " + str(node_indices[quads[k,1]-1])  + " " + str(node_indices[quads[k,2]-1]) + " " + str(node_indices[quads[k,3]-1]) + " 1\n")
            else:
                struct.write(str(node_indices[quads[k,0]-1]) + " " + str(node_indices[quads[k,1]-1])  + " " + str(node_indices[quads[k,2]-1]) + " " + str(node_indices[quads[k,3]-1]) + " 0\n")

    struct.write('\nEnd\n')
    struct.close()
