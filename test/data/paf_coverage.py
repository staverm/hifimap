#!/usr/bin/env python3
import sys

minimap_file = open(sys.argv[1], 'r')
hifimap_file = open(sys.argv[2], 'r')
mismatches = open(sys.argv[3], 'w')

k = 30
if(len(sys.argv) == 5): k = int(sys.argv[4])

minimap = {}
close_mappings = 0

def close(a, b, k = k):
    if abs(int(a) - int(b)) > k: return False
    return True

minimap_lines = minimap_file.readlines()
hifimap_lines = hifimap_file.readlines()

for line in minimap_lines:
   s_line = line.split()
   if len(s_line) < 9: continue
   positions = (s_line[2], s_line[3], s_line[7], s_line[8])
   if s_line[0] not in minimap:
       minimap[s_line[0]] = list()
   minimap[s_line[0]].append(positions)

for line in hifimap_lines:
    if line == False:
        continue
    s_line = line.split()
    if len(s_line) < 9: continue

    if s_line[0] not in minimap: continue
    for start_q, end_q, start_r, end_r in minimap[s_line[0]]:
        if close(start_q, s_line[2]) and close(end_q, s_line[3]) and close(start_r, s_line[7]) and close(end_r, s_line[8]):
            close_mappings += 1
            minimap[s_line[0]].remove((start_q, end_q, start_r, end_r))
            break

mismatches.write("Found in minimap, but not in hifimap: \n")
for read, mappings in minimap.items():
    if len(mappings) > 0:
        mismatches.write(read + '\n')
        for start_q, end_q, start_r, end_r in mappings:
            mismatches.write('\t' + ' ' + start_q + ' ' + end_q + ' ' + start_r + ' ' + end_r + '\n')

print("Minimap mappings: ", len(minimap_lines))
print("Hifimap mappings: ", len(hifimap_lines))
print("Close mappings (k = ", k, ") ", close_mappings)
