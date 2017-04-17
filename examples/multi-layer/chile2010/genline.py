fh = open('line.xyz', 'w')
for i in xrange(0, -61, -1):
    for j in xrange(-120, -59):
        if j in xrange(-100, -90) and i in[-30, -29]:
            txt = '%.1f %.1f %.1f\n' %(j, i, 4.0)
            fh.write(txt)
        else:
            txt = '%.1f %.1f %.1f\n' %(j, i, 0.0)
            fh.write(txt)
fh.close()


fh = open('my_hump.xyz', 'w')

# from numpy import where, array, meshgrid, exp
# x = array(range(-120, -59))
# y = array(range(-60, 0))
# X, Y = meshgrid(x, y)
# ze = -((X+100e0)**2 + (Y+40e0)**2)/10.
# z = where(ze>-10., 40.e0*exp(ze), 0.)
# for i, x1 in enumerate(x):
#     for j, y1 in enumerate(y):
#         txt = '%.1f %.1f %.1f\n' %(x1, y1, ze[i-1, j-1])
#         fh.write(txt)
# fh.close()