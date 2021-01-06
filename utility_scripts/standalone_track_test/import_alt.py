from nucleus_track_3d_2 import *

trk_input = "/mnt/c/Users/hkqur/OneDrive/PHD/pygui/datasets/wtpii/appendcsv.csv"
trk_output= "/mnt/c/Users/hkqur/OneDrive/PHD/pygui/test_saturnia_track/"
frame_no = 200
scan_rad = 22
min_scanR = 19
decay_rate = 1
mit_buff = 30

track_nuc(trk_input,trk_output,frame_no, scan_rad, min_scanR, decay_rate, mit_buff)