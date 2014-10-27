from matplotlib import pyplot as plt
plt.switch_backend("Qt4Agg") 
import numpy as np
from matplotlib_venn import venn2,venn2_circles
from matplotlib_venn import venn3, venn3_circles


figure,axes=plt.subplots(3,2)
#========================= draw GS 3C9, 4B7
v1 = venn2(subsets={'01':220,'10':1268,'11':214},set_labels=('GS_VS_3C9','GS_VS_4B7'),ax=axes[0][0])
v1.get_patch_by_id('01').set_alpha(0.49)
v1.get_patch_by_id('10').set_alpha(0.49)
v1.get_patch_by_id('11').set_alpha(0.49)

c1 = venn2_circles(subsets={'01':220,'10':1268,'11':214},
                  linestyle='dashed',ax=axes[0][0])
first = c1[0]
first.set_edgecolor('red')

second = c1[1]
second.set_edgecolor('green')


#========================== draw GS 91-1C8,91-2A6
v1 = venn2(subsets={'01':1105,'10':519,'11':1646},set_labels=('GS_VS_1C8','GS_VS_2A6'),ax=axes[0][1])
v1.get_patch_by_id('01').set_alpha(0.49)
v1.get_patch_by_id('10').set_alpha(0.49)
v1.get_patch_by_id('11').set_alpha(0.49)

c1 = venn2_circles(subsets={'01':1105,'10':519,'11':1646},
                  linestyle='dashed',ax=axes[0][1])
first = c1[0]
first.set_edgecolor('red')

second = c1[1]
second.set_edgecolor('green')


#=========================== draw GS 3C9-7,8,9
v = venn3(subsets={'001':93,'010':256,'011':93,'100':75,'101':89,'110':72,'111':754}, 
          set_labels = ('GS_VS_3C9-7','GS_VS_3C9-8','GS_VS_3C9-9'),ax=axes[1][0])

v.get_patch_by_id('001').set_alpha(0.49)
v.get_patch_by_id('010').set_alpha(0.49)
v.get_patch_by_id('011').set_alpha(0.49)
v.get_patch_by_id('100').set_alpha(0.49)
v.get_patch_by_id('101').set_alpha(0.49)
v.get_patch_by_id('110').set_alpha(0.49)
v.get_patch_by_id('111').set_alpha(0.49)

label = v.get_label_by_id('011')
label.set_x(label.get_position()[0] + 0.05)

label = v.get_label_by_id('101')
label.set_x(label.get_position()[0] - 0.03)

label = v.get_label_by_id('110')
label.set_x(label.get_position()[0] - 0.1)
label.set_y(label.get_position()[1] + 0.05)

c = venn3_circles(subsets={'001':93,'010':256,'011':93,'100':75,'101':89,'110':72,'111':754},
                  linestyle='dashed',ax=axes[1][0])
first = c[0]
first.set_edgecolor('red')

second = c[1]
second.set_edgecolor('green')

third = c[2]
third.set_edgecolor('blue')





#============================ draw GS 4B7-10,11,12
v = venn3(subsets={'001':86,'010':66,'011':29,'100':46,'101':26,'110':24,'111':199}, 
          set_labels = ('GS_VS_4B7-10','GS_VS_4B7-11','GS_VS_4B7-12'),ax=axes[1][1])

v.get_patch_by_id('001').set_alpha(0.49)
v.get_patch_by_id('010').set_alpha(0.49)
v.get_patch_by_id('011').set_alpha(0.49)
v.get_patch_by_id('100').set_alpha(0.49)
v.get_patch_by_id('101').set_alpha(0.49)
v.get_patch_by_id('110').set_alpha(0.49)
v.get_patch_by_id('111').set_alpha(0.49)

label = v.get_label_by_id('011')
label.set_x(label.get_position()[0] + 0.05)

label = v.get_label_by_id('101')
label.set_x(label.get_position()[0] - 0.03)

label = v.get_label_by_id('110')
label.set_x(label.get_position()[0])
label.set_y(label.get_position()[1])

c = venn3_circles(subsets={'001':86,'010':66,'011':29,'100':46,'101':26,'110':24,'111':199},
                  linestyle='dashed',ax=axes[1][1])
first = c[0]
first.set_edgecolor('red')

second = c[1]
second.set_edgecolor('green')

third = c[2]
third.set_edgecolor('blue')

#=========================== draw GS 91-1C8 13,14,15
v = venn3(subsets={'001':599,'010':213,'011':73,'100':164,'101':112,'110':201,'111':1116}, 
          set_labels = ('GS_VS_1C8-13','GS_VS_1C8-14','GS_VS_1C8-15'),ax=axes[2][0])

v.get_patch_by_id('001').set_alpha(0.49)
v.get_patch_by_id('010').set_alpha(0.49)
v.get_patch_by_id('011').set_alpha(0.49)
v.get_patch_by_id('100').set_alpha(0.49)
v.get_patch_by_id('101').set_alpha(0.49)
v.get_patch_by_id('110').set_alpha(0.49)
v.get_patch_by_id('111').set_alpha(0.49)

label = v.get_label_by_id('011')
label.set_x(label.get_position()[0] + 0.05)

label = v.get_label_by_id('101')
label.set_x(label.get_position()[0] - 0.03)

label = v.get_label_by_id('110')
label.set_x(label.get_position()[0])
label.set_y(label.get_position()[1])

c = venn3_circles(subsets={'001':599,'010':213,'011':73,'100':164,'101':112,'110':201,'111':1116},
                  linestyle='dashed',ax=axes[2][0])
first = c[0]
first.set_edgecolor('red')

second = c[1]
second.set_edgecolor('green')

third = c[2]
third.set_edgecolor('blue')

#=========================== draw GS 91-2A6 16,17,18
v = venn3(subsets={'001':261,'010':119,'011':147,'100':257,'101':88,'110':227,'111':1480}, 
          set_labels = ('GS_VS_2A6-16','GS_VS_2A6-17','GS_VS_2A6-18'),ax=axes[2][1])

v.get_patch_by_id('001').set_alpha(0.49)
v.get_patch_by_id('010').set_alpha(0.49)
v.get_patch_by_id('011').set_alpha(0.49)
v.get_patch_by_id('100').set_alpha(0.49)
v.get_patch_by_id('101').set_alpha(0.49)
v.get_patch_by_id('110').set_alpha(0.49)
v.get_patch_by_id('111').set_alpha(0.49)

label = v.get_label_by_id('011')
label.set_x(label.get_position()[0] + 0.05)

label = v.get_label_by_id('101')
label.set_x(label.get_position()[0] - 0.03)

label = v.get_label_by_id('110')
label.set_x(label.get_position()[0])
label.set_y(label.get_position()[1])

c = venn3_circles(subsets={'001':261,'010':119,'011':147,'100':257,'101':88,'110':227,'111':1480},
                  linestyle='dashed',ax=axes[2][1])
first = c[0]
first.set_edgecolor('red')

second = c[1]
second.set_edgecolor('green')

third = c[2]
third.set_edgecolor('blue')


plt.show()
