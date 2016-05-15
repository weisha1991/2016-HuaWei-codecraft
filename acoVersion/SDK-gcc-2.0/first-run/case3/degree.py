

f = open("topo.csv", "r")
f1 = open("demand.csv", "r")
src_list = []
desc_list = []
demand = []
for line in f1:
    demand = line.strip().split(',')[2].split('|')
print demand
for line in f:
    src_list.append(line.strip().split(',')[1])
    desc_list.append(line.strip().split(',')[2])

cnt=1
for a, b in zip(src_list, desc_list):
    if(a in demand or b in demand):
     	cnt+=1   
print cnt
# vertex_num = max(max(src_list), max(desc_list)) + 1

# out_cnt = [0] * vertex_num
# for node in src_list:
#     out_cnt[node] += 1
# # for i in range(0, len(out_cnt)):
# #     print "%d'out degree=%d" % (i, out_cnt[i])

# in_cnt = [0] * vertex_num
# for node in desc_list:
#     in_cnt[node] += 1
# # for i in range(0, len(in_cnt)):
# #     print "%d'in degree=%d" % (i, in_cnt[i])

# print "*" * 20

# print "out degree 0:"
# for i in range(0, len(in_cnt)):
#     if(out_cnt[i] == 0):
#         print "%d," % (i),

# print "\n"
# print "in degree 0:"
# for i in range(0, len(in_cnt)):
#     if(in_cnt[i] == 0):
#         print "%d," % (i),
# print "\n"
