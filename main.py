import numpy as np
import matplotlib.pyplot as plt
import RIRSimulator.roomsimove_single as RIR
import VEM as vem

rt60 = 0.4 # in seconds
room_dim = [4.2, 3.4, 5.2] # in meters
mic_pos1 = [2, 2, 2] # in  meters
#mic_pos2 = [2, 2, 1] # in  meters
source_pos = [1, 1, 1] # in  meters
sampling_rate = 16000

absorption = RIR.rt60_to_absorption(room_dim, rt60)

# mic_positions = [mic_pos1]
# rir = RIR.do_everything(room_dim, mic_positions, source_pos, rt60)
# np.save('h.npy', rir)
rir = np.load('h.npy')


vem = vem.VEM(rir, 10**(-8), sampling_rate, 200, absorption)

vfe = vem.vfe()
print(vfe)
# plt.plot(rir)
# plt.show()
