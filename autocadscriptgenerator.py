import math

f=open('/home/soenke/Documents/python/script.scp','a');
R=250

for cnt in range(30):
    f.write('3DFACE\r\n')
    xtemp1=round(R*math.cos(math.radians(12*cnt)))
    ytemp1=round(R*math.sin(math.radians(12*cnt)))
    xtemp2=round(R*math.cos(math.radians(12*(cnt+1))))
    ytemp2=round(R*math.sin(math.radians(12*(cnt+1))))
    f.write(str(xtemp1)+','+str(ytemp1)+',-10000\r\n')
    f.write(str(xtemp1)+','+str(ytemp1)+',10000\r\n')
    f.write(str(xtemp2)+','+str(ytemp2)+',10000\r\n')
    f.write(str(xtemp2)+','+str(ytemp2)+',-10000\r\n\r\n')
f.close()
