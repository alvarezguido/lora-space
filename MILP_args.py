import numpy as np
import math


scenario = "wider_scenario_2"
path = "./" + scenario + "/"

#####################
# Positions
#####################

# leo_pos[i,j]   | [time i, x-y-z position]
leo_pos = np.loadtxt(path + "LEO-XYZ-Pos.csv", skiprows=1,
                     delimiter=',', usecols=(1, 2, 3))
# sites_pos[i,j] | [node i, x-y-z position]
sites_pos = np.loadtxt(path + "SITES-XYZ-Pos.csv",
                       skiprows=1, delimiter=',', usecols=(1, 2, 3))

#####################
# Distance
#####################

dist_sat = np.zeros((sites_pos.shape[0], 3, leo_pos.shape[0]))
t = 0
for i in range(leo_pos.shape[0]):
    t += 1
    dist_sat[:, :, i] = leo_pos[i, :] - sites_pos

# distance[i,j] | [node i, distance]
distance = np.zeros((sites_pos.shape[0], leo_pos.shape[0]))
distance[:, :] = (dist_sat[:, 0, :]**2 + dist_sat[:, 1, :]
                  ** 2 + dist_sat[:, 2, :]**2)**(1/2)

#####################
# Rx Power
#####################

freq = 868e6
Lpl = np.zeros((sites_pos.shape[0], leo_pos.shape[0]))
Lpl = 20*np.log10(distance*1000) + 20*np.log10(freq) - \
    147.55  

Ptx = 14
G_device = 0 
G_sat = 12   
Prx = np.zeros((sites_pos.shape[0], leo_pos.shape[0]))
# Prx[i,j] | [node i, prx]
Prx = Ptx + G_sat + G_device - 20 * \
    np.log10(distance*1000) - 20*np.log10(freq) + 147.55
# rssi_nodeid = Prx[nodeid,:]

#####################
# SF
#####################

sf7 = np.array([7, -123, -120, -117.0])
sf8 = np.array([8, -126, -123, -120.0])
sf9 = np.array([9, -129, -126, -123.0])
sf10 = np.array([10, -132, -129, -126.0])
sf11 = np.array([11, -134.53, -131.52, -128.51])
sf12 = np.array([12, -137, -134, -131.0])

spr_fact = np.zeros((sites_pos.shape[0], leo_pos.shape[0]))

for nod_id in range(0, Prx.shape[0]):
    for t_stamp in range(0, Prx.shape[1]):
        if Prx[nod_id, t_stamp] > sf7[1]:
            spr_fact[nod_id, t_stamp] = 7
        elif Prx[nod_id, t_stamp] > sf8[1]:
            spr_fact[nod_id, t_stamp] = 8
        elif Prx[nod_id, t_stamp] > sf9[1]:
            spr_fact[nod_id, t_stamp] = 9
        elif Prx[nod_id, t_stamp] > sf10[1]:
            spr_fact[nod_id, t_stamp] = 10
        elif Prx[nod_id, t_stamp] > sf11[1]:
            spr_fact[nod_id, t_stamp] = 11
        elif Prx[nod_id, t_stamp] > sf12[1]:
            spr_fact[nod_id, t_stamp] = 12
        else:
            spr_fact[nod_id, t_stamp] = 13

#####################
# Airtime
#####################

def airtime(sf):
    pl = 20      # packet len
    cr = 1       # coding rate
    bw = 125     # bandwith
    H = 0        # implicit header disabled (H=0) or not (H=1)
    DE = 0       # low data rate optimization enabled (=1) or not (=0)
    Npream = 8   # number of preamble symbol (12.25  from Utz paper)

    if bw == 125 and sf in [11, 12]:
        DE = 1 # low data rate optimization mandated for BW125 with SF11 and SF12
    if sf == 6:
        H = 1 # can only have implicit header with SF6

    Tsym = (2.0**sf)/bw
    Tpream = (Npream + 4.25)*Tsym
    payloadSymbNB = 8 + max(math.ceil((8.0*pl-4.0*sf+28+16-20*H)/(4.0*(sf-2*DE)))*(cr+4),0)
    Tpayload = payloadSymbNB * Tsym
    return ((Tpream + Tpayload)/1000) # in secs

# for sf in [7,8,9,10,11,12]:
#     print(airtime(sf)*1000)

#####################
# Optimizer
#####################

#proptime = distance[nodeid,:]*(1/c)
max_proptime = 3000 / 299792.458 # in secs

import gurobipy as gp
from gurobipy import GRB
from gurobipy import LinExpr
import os

def solve_milp(nrNodes, packetsToSend, channels, maxBSReceives, sub_channels):

    # TODO: maxBSRx: never hits 16 (3 channels -> 3 concurrent tx)

    # packetsToSend = 2 
    # nodes = [0,1,2,3]
    # beacons = [0,1,2,3]
    # channels = [0,1,2]
    # sfs = [7,8,9,10,11,12]
    # sfm = [ [ 0,  0, 11, 12], # sfm[n][b]
    #         [ 0, 12, 11, 0 ],
    #         [12, 10,  7, 12],
    #         [12, 11, 9, 12]]

    nodes = range(0,nrNodes)
    beacon_dur_secs = 120
    beacons = range(0, int(leo_pos.shape[0]/beacon_dur_secs))
    channels = range(0, channels * sub_channels)
    sfs = [7,8,9,10,11,12]
    sfm = []
    for n,node in enumerate(nodes):
        sfm.append([])
        for b,beacon in enumerate(beacons):
            t_stamp_str = b * beacon_dur_secs
            t_stamp_end = (b + 1) * beacon_dur_secs
            highest_sf = 0
            for t_stamp in range(t_stamp_str, t_stamp_end):
                if spr_fact[n % 1500, t_stamp] > highest_sf and spr_fact[n % 1500, t_stamp] != 13:
                    highest_sf = int(spr_fact[n % 1500, t_stamp])
            sfm[n].append(highest_sf)

    try:
        # Create a new model
        m = gp.Model("lora")

        # Create variables
        var_packet = []
        for n,node in enumerate(nodes):
            var_packet.append([])
            for b,beacon in enumerate(beacons):
                var_packet[n].append([])
                for c,channel in enumerate(channels):
                    var_packet[n][b].append([])
                    for s,sf in enumerate(sfs):
                        var_name = "p_n{}_b{}_c{}_sf{}".format(n,b,channel,sf)
                        var = m.addVar(vtype=GRB.BINARY, name=var_name)
                        var_packet[n][b][c].append(var)
                        # var_packet[n][b][c][s] = var_packet[n][b][c+1][sf-7]

        # Set objective
        obj = LinExpr()
        for n,node in enumerate(nodes):
            for b,beacon in enumerate(beacons):
                for c,channel in enumerate(channels):
                    for s,sf in enumerate(sfs):
                        obj += var_packet[n][b][c][s] * (20 - airtime(sf))
        m.setObjective(obj, GRB.MAXIMIZE)

        # Constraint 1
        for n,node in enumerate(nodes):
            for b,beacon in enumerate(beacons):
                cnt1_name = "c1_n{}_b{}".format(n,b)
                cnt1 = LinExpr()
                for c,channel in enumerate(channels):
                    for s,sf in enumerate(sfs):
                        cnt1 += var_packet[n][b][c][s]
                m.addConstr(cnt1 <= 1, cnt1_name)

        # Constraint 2
        for n,node in enumerate(nodes):
            cnt2_name = "c2_n{}".format(n)
            cnt2 = LinExpr()
            for b,beacon in enumerate(beacons):
                for c,channel in enumerate(channels):
                    for s,sf in enumerate(sfs):
                        cnt2 += var_packet[n][b][c][s]
            m.addConstr(cnt2 <= packetsToSend, cnt2_name)

        # Constraint 3
        for n,node in enumerate(nodes):
            for b,beacon in enumerate(beacons):
                for s,sf in enumerate(sfs):
                    cnt3_name = "c3_n{}_b{}_s{}".format(n,b,sf)
                    cnt3 = LinExpr()
                    for c,channel in enumerate(channels):
                        cnt3 += var_packet[n][b][c][s]
                    m.addConstr(cnt3 <= 1, cnt3_name)

        # Constraint 4
        for b,beacon in enumerate(beacons):
            for c,channel in enumerate(channels):
                cnt4_name = "c4_b{}_c{}".format(b,c)
                cnt4 = LinExpr()
                for n,node in enumerate(nodes):
                    for s,sf in enumerate(sfs):
                        cnt4 += var_packet[n][b][c][s] * (airtime(sf) + max_proptime)
                m.addConstr(cnt4 <= beacon_dur_secs, cnt4_name)

        # Constraint 5
        M = 100
        for n,node in enumerate(nodes):
            for b,beacon in enumerate(beacons):
                for c,channel in enumerate(channels):
                    for s,sf in enumerate(sfs):
                        cnt5_name = "c5_n{}_b{}_c{}_sf{}".format(n,b,c,sf)
                        cnt5 = LinExpr()
                        if sfm[n][b] == 0:
                            cnt5 += var_packet[n][b][c][s]
                        else:
                            cnt5 += sfm[n][b] - sf
                        m.addConstr(cnt5 <= M * (1 - var_packet[n][b][c][s]), cnt5_name)

        # Optimize model
        m.Params.TimeLimit = 60 # 1 min
        m.optimize()
        # print('Obj: %g' % m.objVal)
        
        # Save log
        logs = []
        time = 0
        for b,beacon in enumerate(beacons):
            beacon_pcount = 0
            beacon_airtime = 0
            time = 120 * b
            logs.append("{:3.3f},B,[]".format(time + 2))
            for c,channel in enumerate(channels):
                for n,node in enumerate(nodes):
                    for s,sf in enumerate(sfs):
                        if var_packet[n][b][c][s].x != 0:
                            #print('%s %g' % (var_packet[n][b][c][s].varName, var_packet[n][b][c][s].x))   
                            beacon_pcount += 1
                            beacon_airtime += airtime(sf) + max_proptime
                            time += airtime(sf) + max_proptime
                            if time > 1200:
                                time = 1200
                            logs.append("{:3.3f},{},{:3.3f},{:3.3f},{},PE".format(time, n, distance[n % 1500,math.ceil(time)], 0, sf))
            # print('beacon {} - p{} a{:.2f}'.format(b, beacon_pcount, beacon_airtime))
    
    except gp.GurobiError as e:
        print('Error code ' + str(e.errno) + ': ' + str(e))
    except AttributeError:
        print('Encountered an attribute error')

    return logs

#####################
# Solve
#####################

# LoRa-L
# name = 'LO'
# node_num = [5, 10, 25, 50, 100, 250, 500, 1000, 1400, 2000, 3000, 4000, 5000, 6000]
# pack_num = [1, 2, 3, 5, 8, 9, 10]
# maxBSReceives = 16
# sub_channels = 1

# LoRa-E
name = 'EO'
node_num = [50, 500, 750, 1000, 2000, 4000, 6000, 8000, 10000, 12000, 14000, 16000, 18000, 20000]
pack_num = [3]
maxBSReceives = 500
sub_channels = 280

RANDOM_SEED = 0
channels = [1, 3]

for nrNodes in node_num:
    for packetsToSend in pack_num:
        for channel in channels:

            logs = solve_milp(nrNodes, packetsToSend, channel, maxBSReceives, sub_channels)

            folder = "./" + name + '_' + str(channel) + 'CH_s' + str(RANDOM_SEED) + '_p' + str(packetsToSend)
            if not os.path.exists(folder):
                os.makedirs(folder)

            fname = folder + "/" + name + "_p" + str(packetsToSend) + "_" + str(nrNodes) + "_" + str(channel)+"CH_" + str(maxBSReceives) + ".csv"

            with open(fname,"w") as myfile:
                myfile.write("\n".join(logs))
            myfile.close()
            print('SOLVED:', fname)

