#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Oct  2 11:08:09 2020

@author: alvarezguido
GITHUB: https://github.com/alvarezguido
"""

"""
SYNOPSIS
----
----
-----

"""

import simpy
import random
import numpy as np
import math
#import sys
#import re
import matplotlib.pyplot as plt
#import os
#import operator
from mpl_toolkits import mplot3d
from mpl_toolkits.mplot3d import Axes3D
#import PIL
import random
import re
import os
import datetime
import sys

name = "LTb"
mode_debbug = 0

beacon_rec = 0
max_rec = 15

if not mode_debbug:
    null = open(os.devnull, 'w')
    old_stdout = sys.stdout
    sys.stdout = null

####WE START BY USING SF=12 ADN BW=125 AND CR=1, FOR ALL NODES AND ALL TRANSMISIONS######
if mode_debbug:
    RANDOM_SEED = 5
    chan = 1
    packetlen = 20
    total_data = 60
    beacon_time = 120
    maxBSReceives = 16
    multi_nodes = [10]
    p_skip_param = 4000
else:
    RANDOM_SEED = int(sys.argv[1])
    chan = int(sys.argv[2])
    packetlen = int(sys.argv[3])   ##NODES SEND PACKETS OF JUST 20 Bytes
    total_data = int(sys.argv[4]) ##TOTAL DATA ON BUFFER, FOR EACH NODE (IT'S THE BUFFER O DATA BEFORE START SENDING)
    beacon_time = int(sys.argv[5]) ###SAT SENDS BEACON EVERY CERTAIN TIME
    maxBSReceives = int(sys.argv[6]) ##MAX NUMBER OF PACKETS THAT BS (ie SATELLITE) CAN RECEIVE AT SAME TIME
    
    multi_nodes = [int(sys.argv[7]), int(sys.argv[8]) ,int(sys.argv[9]), int(sys.argv[10]),int(sys.argv[11]),int(sys.argv[12]),int(sys.argv[13]),int(sys.argv[14]),int(sys.argv[15]),int(sys.argv[16]),int(sys.argv[17]),int(sys.argv[18]),int(sys.argv[19]),int(sys.argv[20])]
    p_skip_param = int(sys.argv[21])

random.seed(RANDOM_SEED) #RANDOM SEED IS FOR GENERATE ALWAYS THE SAME RANDOM NUMBERS (ie SAME RESULTS OF SIMULATION)
nodesToSend = []
packetsToSend = math.ceil(total_data/packetlen)

###GLOBAL PARAMS ####
bsId = 1 ##ID OF BASE STATION (NOT USED)
channel = [0,1,2] ##NOT USED BY NOW

avgSendTime = 3  ## NOT USED! --> A NODE SENDS A PACKET EVERY X SECS
back_off = beacon_time * 0.95 ###BACK OFF TIME FOR SEND A PACKET
packetsAtBS = [] ##USED FOR CHEK IF THERE ARE ALREADY PACKETS ON THE SATELLITE
c = 299792.458 ###SPEED LIGHT [km/s]
Ptx = 14
G_device = 0; ##ANTENNA GAIN FOR AN END-DEVICE
G_sat = 12;   ##ANTENNA GAIN FOR SATELLITE
nodes = [] ###EACH NODE WILL BE APPENDED TO THIS VARIABLE
freq =868e6 ##USED FOR PATH LOSS CALCULATION
frequency = [868100000, 868300000, 868500000] ##FROM LORAWAN REGIONAL PARAMETERS EU863-870 / EU868

nrLost = 0 ### TOTAL OF LOST PACKETS DUE Lpl
nrCollisions = 0 ##TOTAL OF COLLIDED PACKETS
nrProcessed = 0 ##TOTAL OF PROCESSED PACKETS
nrReceived = 0 ###TOTAL OF RECEIVED PACKETS


##ARRAY WITH MEASURED VALUES FOR SENSIBILITY, NEW VALUES
##THE FOLLOWING VALUES CORRESPOND TO:
#   - FIRST ELEMENT: IT'S THE SF (NOT USABLE)
#   - SECOND ELEMENT: SENSIBILITY FOR 125KHZ BW
#   - THIRD ELEMENT: SENSIBILITY FOR 250KHZ BW
#   - FOURTH ELEMENT: SENSIBILITY FOR 500KHZ BW
# NOTICE THAT SENSIBILITY DECREASE ALONG BW INCREASES, ALSO WITH LOWER SF
# THIS VALUES RESPONDS TO:
# wf = -174 + 10 log(BW) +NF +SNRf
sf7 = np.array([7,-123,-120,-117.0])
sf8 = np.array([8,-126,-123,-120.0])
sf9 = np.array([9,-129,-126,-123.0])
sf10 = np.array([10,-132,-129,-126.0])
sf11 = np.array([11,-134.53,-131.52,-128.51])
sf12 = np.array([12,-137,-134,-131.0])

sensi = np.array([sf7,sf8,sf9,sf10,sf11,sf12])

path = "./wider_scenario_2/"

### -137dB IS THE MINIMUN TOLERABLE SENSIBILITY, FOR SF=12 AND BW=125KHz ###

leo_pos=np.loadtxt( path + "LEO-XYZ-Pos.csv",skiprows=1,delimiter=',',usecols=(1,2,3))
## WHERE:
    ## leo_pos[i,j]:
        ## i --> the step time in sat pass
        ## j --> 0 for x-position, 1 for y-position, 2 for z-position

sites_pos = np.loadtxt( path + "SITES-XYZ-Pos.csv",skiprows=1,delimiter=',',usecols=(1,2,3))
## WHERE:
    ## sites_pos[i,j]:
        ## i --> the node i
        ## j --> 0 for x-position, 1 for y-position, 2 for z-position


dist_sat = np.zeros((sites_pos.shape[0],3,leo_pos.shape[0]))
t = 0
for i in range(leo_pos.shape[0]):
    t+=1
    dist_sat [:,:,i] = leo_pos[i,:] - sites_pos
## WHERE:
    ## dist_sat[i,j,k]:
        ## i --> the node i
        ## j --> 0 for x-position, 1 for y-position, 2 for z-position
        ## k --> the step time in sat pass
    
#### FOR COMPUTE DISTANCE MAGNITUDE (ABS) FROM END-DEVICE TO SAT PASSING BY ####
distance = np.zeros((sites_pos.shape[0],leo_pos.shape[0]))
distance[:,:] = (dist_sat[:,0,:]**2 + dist_sat[:,1,:]**2 + dist_sat[:,2,:]**2)**(1/2)
## WHERE:
    ## distance[i,j]:
        ## i --> the node i
        ## j --> the step time in sat pass


##MATRIX FOR LINK BUDGET Lpl ###
Lpl = np.zeros((sites_pos.shape[0],leo_pos.shape[0])) 
Lpl = 20*np.log10(distance*1000) + 20*np.log10(freq) - 147.55 #DISTANCE MUST BE IN METERS
## WHERE:
    ## Lpl[i,j]:
        ## i --> the node i
        ## j --> the step time in sat pass 

##MATRIX FOR LINK BUDGET, USING Prx ###
Prx = np.zeros((sites_pos.shape[0],leo_pos.shape[0])) 
Prx = Ptx + G_sat + G_device -20*np.log10(distance*1000) - 20*np.log10(freq) + 147.55 #DISTANCE IS CONVERTED TO METERS
## WHERE:
    ## Prx[i,j]:
        ## i --> the node i
        ## j --> the step time in sat pass 
distance = np.concatenate((distance,distance,distance,distance,distance,distance,distance,distance,distance,distance,distance,distance,distance,distance,distance,distance,distance,distance,distance,distance))
Lpl = np.concatenate((Lpl,Lpl,Lpl,Lpl,Lpl,Lpl,Lpl,Lpl,Lpl,Lpl,Lpl,Lpl,Lpl,Lpl,Lpl,Lpl,Lpl,Lpl,Lpl,Lpl))
Prx = np.concatenate((Prx,Prx,Prx,Prx,Prx,Prx,Prx,Prx,Prx,Prx,Prx,Prx,Prx,Prx,Prx,Prx,Prx,Prx,Prx,Prx))

elev = np.degrees(np.arcsin(599/distance))

IS7 = np.array([1,-8,-9,-9,-9,-9])
IS8 = np.array([-11,1,-11,-12,-13,-13])
IS9 = np.array([-15,-13,1,-13,-14,-15])
IS10 = np.array([-19,-18,-17,1,-17,-18])
IS11 = np.array([-22,-22,-21,-20,1,-20])
IS12 = np.array([-25,-25,-25,-24,-23,1])

#THIS IS THE MATRIX OF CROSS INTERFERENCE BETWEEN SF
IsoThresholds = np.array([IS7,IS8,IS9,IS10,IS11,IS12])
Collmap = [[0 for i in range(0,6)] for j in range(0,6)]

def simulate_scenario (nrNodes):
    env = simpy.Environment()
    
    # check only the capture between the same spreading factor
    def powerCollision_1(p1, p2):
        #powerThreshold = 6
        print ("pwr: node {} with rssi {} dBm and node {} with rssi {} dBm; diff {:3.2f} dBm".format(p1,p1.rssi[math.ceil(env.now)], p2.nodeid,p2.rssi[math.ceil(env.now)], p1.rssi[math.ceil(env.now)] - p2.rssi[math.ceil(env.now)]))
        if p1.sf == p2.sf:
           if abs(p1.rssi[math.ceil(env.now)] - p2.rssi[math.ceil(env.now)]) < IsoThresholds[p1.sf-7][p2.sf-7]:
                print ("collision pwr both node {} and node {}".format(p1.nodeid, p2.nodeid))
                # packets are too close to each other, both collide
                # return both pack ets as casualties
                return (p1, p2)
           elif p1.rssi[math.ceil(env.now)] - p2.rssi[math.ceil(env.now)] < IsoThresholds[p1.sf-7][p2.sf-7]:
                # p2 overpowered p1, return p1 as casualty
                print ("collision pwr node {} overpowered node {}".format(p2.nodeid, p1.nodeid))
                return (p1,)
           print ("p1 wins, p2 lost")
           # p2 was the weaker packet, return it as a casualty
           return (p2,)
        else:
           return ()
    
    # check the capture effect and checking the effect of pesudo-orthognal SFs
    def powerCollision_2(p1, p2):
        #powerThreshold = 6
        global Collmap
        print ("SF: node {} has {} ; node {} has {}".format(p1.nodeid,p1.sf, p2.nodeid, p2.sf))
        print ("pwr: node {} with rssi {} dBm and node {} with rssi {} dBm; diff {:3.2f} dBm".format(p1,p1.rssi[math.ceil(env.now)], p2.nodeid,p2.rssi[math.ceil(env.now)], p1.rssi[math.ceil(env.now)] - p2.rssi[math.ceil(env.now)]))
        if p1.sf == p2.sf:
           if abs(p1.rssi[math.ceil(env.now)] - p2.rssi[math.ceil(env.now)]) < IsoThresholds[p1.sf-7][p2.sf-7]:
               print ("collision pwr both node {} and node {}".format(p1.nodeid, p2.nodeid))
               Collmap[p1.sf-7][p2.sf-7] += 1
               Collmap[p2.sf-7][p1.sf-7] += 1
               # packets are too close to each other, both collide
               # return both packets as casualties
               return (p1, p2)
           elif p1.rssi[math.ceil(env.now)] - p2.rssi[math.ceil(env.now)] < IsoThresholds[p1.sf-7][p2.sf-7]:
               # p2 overpowered p1, return p1 as casualty
               print ("collision pwr node {} overpowered node {}".format(p2.nodeid, p1.nodeid))
               print ("capture - p2 wins, p1 lost")
               Collmap[p1.sf-7][p2.sf-7] += 1
               return (p1,)
           print ("capture - p1 wins, p2 lost")
           # p2 was the weaker packet, return it as a casualty
           Collmap[p2.sf-7][p1.sf-7] += 1
           return (p2,)
        else:
           if p1.rssi[math.ceil(env.now)]-p2.rssi[math.ceil(env.now)] > IsoThresholds[p1.sf-7][p2.sf-7]:
              print ("P1 is OK")
              if p2.rssi[math.ceil(env.now)]-p1.rssi[math.ceil(env.now)] > IsoThresholds[p2.sf-7][p1.sf-7]:
                  print ("p2 is OK")
                  return ()
              else:
                  print ("p2 is lost")
                  Collmap[p2.sf-7][p1.sf-7] += 1
                  return (p2,)
           else:
               print ("p1 is lost")
               Collmap[p1.sf-7][p2.sf-7] += 1
               if p2.rssi[math.ceil(env.now)]-p1.rssi[math.ceil(env.now)] > IsoThresholds[p2.sf-7][p1.sf-7]:
                   print ("p2 is OK")
                   return (p1,)
               else:
                   print( "p2 is lost")
                   Collmap[p2.sf-7][p1.sf-7] += 1
                   return (p1,p2)
        
       
    def timingCollision(p1, p2):
        # assuming p1 is the freshly arrived packet and this is the last check
        # we've already determined that p1 is a weak packet, so the only
        # way we can win is by being late enough (only the first n - 5 preamble symbols overlap)
    
        # assuming 8 preamble symbols
        Npream = 8
    
        # we can lose at most (Npream - 5) * Tsym of our preamble
        Tpreamb = 2**p1.sf/(1.0*p1.bw) * (Npream - 5)
    
        # check whether p2 ends in p1's critical section
        p2_end = p2.addTime + p2.rectime
        p1_cs = env.now + (Tpreamb/1000.0)  # to sec
        print ("collision timing node {} ({},{},{}) node {} ({},{})".format(
            p1.nodeid, env.now - env.now, p1_cs - env.now, p1.rectime,
            p2.nodeid, p2.addTime - env.now, p2_end - env.now
        ))
        if p1_cs < p2_end:
            # p1 collided with p2 and lost
            print ("not late enough")
            return True
        print ("saved by the preamble")
        return False
    
    def checkcollision(packet):
        col = 0 # flag needed since there might be several collisions for packet
        processing = 0
        #print ("MAX RECEIVE IS: ", maxBSReceives)
        for i in range(0,len(packetsAtBS)):
            if packetsAtBS[i].packet.processed == 1:
                processing = processing + 1
        if (processing > maxBSReceives):
            print ("{:3.5f} || Too much packets on Base Sattion.. Packet will be lost!", len(packetsAtBS))
            packet.processed = 0
        else:
            packet.processed = 1
    
        if packetsAtBS:
            print ("{:3.5f} || >> FOUND overlap... node {} (sf:{} bw:{} freq:{}) others: {}".format(env.now,packet.nodeid, packet.sf, packet.bw,packet.freq,len(packetsAtBS)))
            for other in packetsAtBS:
                if other.nodeid != packet.nodeid:
                   print ("{:3.5f} || >> node {} overlapped with node {} (sf:{} bw:{} freq:{}). Let's check Freq...".format(env.now,packet.nodeid, other.nodeid, other.packet.sf, other.packet.bw,other.packet.freq))
                   # simple collision
                   #if frequencyCollision(packet, other.packet) and sfCollision(packet, other.packet):
                   if frequencyCollision(packet, other.packet) and timingCollision(packet, other.packet):
                       c = powerCollision_2(packet, other.packet)
                       for p in c:
                          p.collided = 1
                          if p == packet:
                             col = 1         
            return col
        return 0
    
    
    ###frequencyCollision, CONDITIONS###
    
    ##|f1-f2| <= 120 kHz if f1 or f2 has bw 500
    ##|f1-f2| <= 60 kHz if f1 or f2 has bw 250
    ##|f1-f2| <= 30 kHz if f1 or f2 has bw 125
    def frequencyCollision(p1,p2):
        if (abs(p1.freq-p2.freq)<=120 and (p1.bw==500 or p2.freq==500)):
            print ("{:3.5f} || >> freq coll on node {} and node {}.. Let's check SF...".format(env.now,p1.nodeid, p2.nodeid))
            return True
        elif (abs(p1.freq-p2.freq)<=60 and (p1.bw==250 or p2.freq==250)):
            print ("{:3.5f} || >> freq coll on node {} and node {}.. Let's check SF...".format(env.now,p1.nodeid, p2.nodeid))
            return True
        else:
            if (abs(p1.freq-p2.freq)<=30):
                print( "{:3.5f} || >> Freq coll on node {} and node {}.. Let's check SF...".format(env.now,p1.nodeid, p2.nodeid))
                return True
            #else:
        print ("{:3.5f} || >> No frequency collision..".format(env.now))
        return False
    
    #FOLLOWING FUNCTION NOT USED
    def channelCollision(p1,p2):
        if (p1.ch == p2.ch):
            print ("{:3.5f} || >> channel coll for ch {} on node {} and ch {} on node {}.. Let's check SF...".format(env.now,p1.ch,p1.nodeid,p2.ch,p2.nodeid))
            return True
        else:
            print ("{:3.5f} || >> No channel collision..".format(env.now))
            return False
    
    def sfCollision(p1, p2):
        if p1.sf == p2.sf:
            print ("{:3.5f} || >> COLLISION! SF coll on node {} and node {} (ie same SF)...".format(env.now,p1.nodeid, p2.nodeid))
            # p2 may have been lost too, will be marked by other checks
            return True
        print ("{:3.5f} || >> No SF Collision!".format(env.now))
        return False
    
     
    def powerCollision(p1, p2):
        powerThreshold = 6 # dB
        print ("{:3.5f} || power: node {} {:3.2f} dBm, node {} {:3.2f}; diff is {}dBm".format(env.now,p1.nodeid,p1.rssi[math.ceil(env.now)],p2.nodeid, p2.rssi[math.ceil(env.now)], round(p1.rssi[math.ceil(env.now)] - p2.rssi[math.ceil(env.now)],2)))
        #print ("pwr: node {0.nodeid} {0.rssi:3.2f} dBm node {1.nodeid} {1.rssi:3.2f} dBm; diff {2:3.2f} dBm".format(p1, p2, round(p1.rssi[math.ceil(env.now)] - p2.rssi[math.ceil(env.now)],2)))
        if abs(p1.rssi[math.ceil(env.now)] - p2.rssi[math.ceil(env.now)]) < powerThreshold:
            print( "{:3.5f} || Collision power both node {} and node {}".format(env.now,p1.nodeid, p2.nodeid))
            # packets are too close to each other, both collide
            # return both packets as casualties
            return (p1, p2)
        elif p1.rssi[math.ceil(env.now)] - p2.rssi[math.ceil(env.now)] < powerThreshold:
            # p2 overpowered p1, return p1 as casualty
            print ("{:3.5f} || Collision pwr node {} has overpowered node {}".format(env.now,p2.nodeid, p1.nodeid))
            return (p1,)
        print ("{:3.5f} || p1 wins, p2 lost".format(env.now))
        # p2 was the weaker packet, return it as a casualty
        return (p2,)
    
    class myNode():
        def __init__(self, nodeid, bs, avgSendTime, packetlen, total_data):
            global channel
            self.nodeid = nodeid
            self.avgSendTime = avgSendTime
            self.bs = bs
            self.dist = distance[nodeid,:]
            self.elev = elev[nodeid,:]
            self.mindist = np.amin(distance[nodeid,:])
            self.mindist_pos = int(np.where(distance[nodeid,:] == np.amin(distance[nodeid,:]))[0])
            self.buffer = total_data
            self.packetlen = packetlen
            #self.ch = int(random.choice(channel)) 
            self.packet = myPacket(self.nodeid, packetlen, self.dist)
            self.sent = 0 #INITIAL SENT PACKETS
            self.totalLost = 0 #INITIAL TOTAL LOST FOR PARTICULAR NODE
            self.totalColl = 0
            self.totalRec = 0
            self.totalProc = 0
            
    
    class myPacket():
        def __init__(self, nodeid, packetlen, dist):
            #global experiment
            global Ptx
            global Prx
            #global gamma
            #global d0
            #global var
            global Lpl
            #global freq
            #global GL
            global c
            global distance
            global channel
            global frequency
            #SF = [7,8,9,10,11,12]
    
            self.nodeid = nodeid
            self.txpow = Ptx
            #self.sf = random.choice(SF)
            self.sf = 12
            self.cr = 1 ##CODING RATE
            self.bw = 125    
            # transmission range, needs update XXX
            self.transRange = 150
            self.pl = packetlen
            self.symTime = (2.0**self.sf)/self.bw
            self.arriveTime = 0
            self.rssi = Prx[nodeid,:]
            self.freq = int(random.choice(frequency)) 
            self.rectime = airtime(self.sf,self.cr,self.pl,self.bw) ##RECTIME IS THE RECEPTION TIME (ie AIRTIME)
            self.proptime = distance[nodeid,:]*(1/c)
            self.collided = 0
            self.processed = 0
            self.lost = bool
    
    
    def airtime(sf,cr,pl,bw):
        H = 0        # implicit header disabled (H=0) or not (H=1)
        DE = 0       # low data rate optimization enabled (=1) or not (=0)
        Npream = 8   # number of preamble symbol (12.25  from Utz paper)
    
        if bw == 125 and sf in [11, 12]:
            # low data rate optimization mandated for BW125 with SF11 and SF12
            DE = 1
        if sf == 6:
            # can only have implicit header with SF6
            H = 1
    
        Tsym = (2.0**sf)/bw
        Tpream = (Npream + 4.25)*Tsym
        #print ("PARAMS FOR TRANSMISION: sf", sf, " cr", cr, "pl", pl, "bw", bw)
        payloadSymbNB = 8 + max(math.ceil((8.0*pl-4.0*sf+28+16-20*H)/(4.0*(sf-2*DE)))*(cr+4),0)
        Tpayload = payloadSymbNB * Tsym
        return ((Tpream + Tpayload)/1000) ##IN SECS
    
    def selectSF (env, node):
        global sf7,sf8,sf9,sf10,sf11,sf12 
        rssi = node.packet.rssi[math.ceil(env.now)]
        #print ("{:3.5f} || RSSI for node {} is {} dB...".format(env.now,node.nodeid,rssi))
        if rssi > sf7[1]:
            #print ("----Select SF7")
            node.packet.sf = 7
        elif rssi > sf8[1]:
            #print ("----Select SF8")
            node.packet.sf = 8
        elif rssi > sf9[1]:
            #print ("----Select SF9")
            node.packet.sf = 9
        elif rssi > sf10[1]:
            #print ("----Select SF10")
            node.packet.sf = 10
        elif rssi > sf11[1]:
            #print ("----Select SF11")
            node.packet.sf = 11
        else:
            #print ("----Select S12")
            node.packet.sf = 12
    
    def transmit(env,node):
        #while nodes[node.nodeid].buffer > 0.0:
        global wait_min
        global wait_max
        global back_off
        global beacon_time
        global logs
        global nodesToSend
        global beacon_rec
        global max_rec
        while node.buffer > 0.0:
            node.packet.sf = 12
            yield env.timeout(node.packet.rectime + float(node.packet.proptime[math.ceil(env.now)])) ##GIVE TIME TO RECEIVE BEACON
                          
            if node in packetsAtBS:
                print ("{:3.5f} || ERROR: packet is already in...".format(env.now))
            else:
                sensibility = sensi[node.packet.sf - 7, [125,250,500].index(node.packet.bw) + 1]
                if node.packet.rssi[math.ceil(env.now)] < sensibility: #HERE WE ARE CONSIDERING RSSI AT TIME ENV.NOW
                    print ("{:3.5f} || Node {}: Can not reach beacon due Lpl".format(env.now,node.nodeid))
                    wait =0 ##LETS WAIT FOR NEXT BEACON
                    node.packet.lost = False
                    trySend = False
    
                else:
                    nodesToSend.append(node.nodeid)
                    beacon_rec += 1
                    wait = random.uniform(1,back_off - node.packet.rectime - float(node.packet.proptime[math.ceil(env.now)])) ##TRIGGER BACK-OFF TIME
                    yield env.timeout(wait)
                    print ("{:3.5f} || Node {} begins to transmit a packet".format(env.now,node.nodeid))
                    selectSF(env,node) ##CHOOSE SF
                    #trySend = True
                    #node.sent = node.sent + 1
                    #node.buffer = node.buffer - node.packetlen
                    if node in packetsAtBS:
                        print ("{} || ERROR: packet is already in...".format(env.now))
                    else:
                        sensibility = sensi[node.packet.sf - 7, [125,250,500].index(node.packet.bw) + 1]
                        if node.packet.rssi[math.ceil(env.now)] < sensibility: #HERE WE ARE CONSIDERING RSSI AT TIME ENV.NOW
                            print ("{:3.5f} || Node {}: The Packet will be Lost due Lpl".format(env.now,node.nodeid))
                            node.packet.lost = True ## LOST ONLY CONSIDERING Lpl
                            trySend = False
                        else:
                            
                            p_skip = 2/(1+math.exp(-beacon_rec/p_skip_param))-1
                            this_p = random.uniform(0,1)
                            #sys.stdout = old_stdout
                            #print('p_skip_param:', p_skip_param, ' beacon_rec:', beacon_rec, ' p_skip:', p_skip, ' this_p:', this_p)
                            if this_p < p_skip:
                                # beacon_rec +=1
                                # nrNodes p_skip_param
                                trySend = 0
                                print ("***********No send!")
                                #sys.stdout = null
                            else:
                                #print ("***********Send!")
                                #sys.stdout = null
                                node.packet.lost = False ## LOST ONLY CONSIDERING Lpl
                                trySend = True
                                node.sent = node.sent + 1
                                node.buffer = node.buffer - node.packetlen
                                print ("{:3.5f} || Prx for node {} is {:3.2f} dB".format(env.now, node.nodeid, node.packet.rssi[math.ceil(env.now)]))
                                #print ("Prx for node",node.nodeid, "is: ",node.packet.rssi[math.ceil(env.now)],"at time",env.now)
                                print ("{:3.5f} || Let's try if there are collisions...".format(env.now))
                                if (checkcollision(node.packet)==1):
                                    node.packet.collided = 1
                                else:
                                    node.packet.collided = 0
                                    print ("{:3.5f} || ...No Collision by now!".format(env.now))
                                packetsAtBS.append(node)
                                node.packet.addTime = env.now
                                yield env.timeout(node.packet.rectime)
                        
            
            if trySend == 1:
                if node.packet.lost:
                    logs.append("{:3.3f},{},{:3.3f},{:3.3f},{},PL".format(env.now,node.nodeid,node.dist[math.ceil(env.now)],node.elev[math.ceil(env.now)],node.packet.sf))
                elif node.packet.collided:
                    logs.append("{:3.3f},{},{:3.3f},{:3.3f},{},PC".format(env.now,node.nodeid,node.dist[math.ceil(env.now)],node.elev[math.ceil(env.now)],node.packet.sf))
                elif node.packet.processed == 0:
                    logs.append("{:3.3f},{},{:3.3f},{:3.3f},{},NP".format(env.now,node.nodeid,node.dist[math.ceil(env.now)],node.elev[math.ceil(env.now)],node.packet.sf))
                else:
                    logs.append("{:3.3f},{},{:3.3f},{:3.3f},{},PE".format(env.now,node.nodeid,node.dist[math.ceil(env.now)],node.elev[math.ceil(env.now)],node.packet.sf))
            
            # complete packet has been received by base station
            # Let's remove from Base Station
            if (node in packetsAtBS):
                packetsAtBS.remove(node)
                # reset the packet
            node.packet.collided = 0
            node.packet.processed = 0
            node.packet.lost = False
            node.packet.sf = 12
            
            #yield env.timeout(beacon_time-wait-node.packet.rectime)
            if trySend:
                yield env.timeout(beacon_time-wait-2*node.packet.rectime)
            else:
                yield env.timeout(beacon_time-wait-node.packet.rectime)
        
                                         
    def beacon (env):
        global beacon_time
        global nodesToSend
        global logs
        global beacon_rec
        i = 0
        while True:
            if i == 0:
                yield env.timeout(0)           
            else:
                yield env.timeout(beacon_time-2)
            i=i+1
            print ("{:3.5f} || ***A new beacon has been sended from Satellite***".format(env.now))
            beacon_rec =0
            yield env.timeout(2)
            logs.append("{:3.3f},B,{}".format(env.now,nodesToSend))
            nodesToSend = []    
        
    env.process(beacon(env)) ##BEACON SENDER
    
    ### THIS IS GOING TO CREATE NODES AND DO TRAMSMISIONS. IS THE MAIN PROGRAM ###
    for i in range(nrNodes):
        node = myNode(i,bsId, avgSendTime, packetlen, total_data)
        nodes.append(node)
        env.process(transmit(env,node))
        
    env.run(until=600*2)
    
    sent = sum(n.sent for n in nodes)
    
    return ([sent,nrCollisions,nrLost,nrProcessed,nrReceived],logs)


#############################################################
if chan == 1:
    ###SCENARIO 1 CHANNEL###
    frequency = [868100000] #1 CH
    
    nodes = [] ###EACH NODE WILL BE APPENDED TO THIS VARIABLE
    nrLost = 0 ### TOTAL OF LOST PACKETS DUE Lpl
    nrCollisions = 0 ##TOTAL OF COLLIDED PACKETS
    nrProcessed = 0 ##TOTAL OF PROCESSED PACKETS
    nrReceived = 0 ###TOTAL OF RECEIVED PACKETS
    
    i =0
    scenario_1ch = np.zeros((len(multi_nodes),5))
    results = []
    ## WHERE:
        ## scenario_1ch[i,j]:
            ## i --> the node i
            ## j --> [sent, nrCollisions, nrLost, nrProcessed, nrReceived]
    
    for nrNodes in multi_nodes:
        print ("\n\n***NEW SCENARIO BEGINS***\n")
        logs = []
        results,logs = simulate_scenario(nrNodes)
        scenario_1ch[i,:] = results
        folder = name+'_1CH_s'+str(RANDOM_SEED)+'_p'+str(packetsToSend)
        if not os.path.exists(folder):
            os.makedirs(folder)
        fname = "./"+folder+"/" + str(name+"_"+str(nrNodes)+"_1CH_"+str(maxBSReceives)+"_s"+str(RANDOM_SEED)+"_p"+str(packetsToSend)) + ".csv"
        with open(fname,"w") as myfile:
            myfile.write("\n".join(logs))
        myfile.close()
        i=i+1
        if not mode_debbug:
            nodes = [] ###EACH NODE WILL BE APPENDED TO THIS VARIABLE
        nrLost = 0 ### TOTAL OF LOST PACKETS DUE Lpl
        nrCollisions = 0 ##TOTAL OF COLLIDED PACKETS
        nrProcessed = 0 ##TOTAL OF PROCESSED PACKETS
        nrReceived = 0 ###TOTAL OF RECEIVED PACKETS


#############################################################
if chan == 3:
    ###SCENARIO 3 CHANNELS###
    frequency = [868100000, 868300000, 868500000] ##FROM LORAWAN REGIONAL PARAMETERS EU863-870 / EU868
    
    nodes = [] ###EACH NODE WILL BE APPENDED TO THIS VARIABLE
    nrLost = 0 ### TOTAL OF LOST PACKETS DUE Lpl
    nrCollisions = 0 ##TOTAL OF COLLIDED PACKETS
    nrProcessed = 0 ##TOTAL OF PROCESSED PACKETS
    nrReceived = 0 ###TOTAL OF RECEIVED PACKETS
    
    i =0
    scenario_3ch = np.zeros((len(multi_nodes),5))
    results = []
    ## WHERE:
        ## scenario_3ch[i,j]:
            ## i --> the node i
            ## j --> [sent, nrCollisions, nrLost, nrProcessed, nrReceived]
    for nrNodes in multi_nodes:
        print ("\n\n***NEW SCENARIO BEGINS***\n")
        logs = []
        results,logs = simulate_scenario(nrNodes)
        scenario_3ch[i,:] = results
        folder = name+'_3CH_s'+str(RANDOM_SEED)+'_p'+str(packetsToSend)
        if not os.path.exists(folder):
            os.makedirs(folder)
        fname = "./"+folder+"/" + str(name+"_"+str(nrNodes)+"_3CH_"+str(maxBSReceives)+"_s"+str(RANDOM_SEED)+"_p"+str(packetsToSend)) + ".csv"
        with open(fname,"w") as myfile:
            myfile.write("\n".join(logs))
        myfile.close()
        i=i+1
        if not mode_debbug:
            nodes = [] ###EACH NODE WILL BE APPENDED TO THIS VARIABLE
        nrLost = 0 ### TOTAL OF LOST PACKETS DUE Lpl
        nrCollisions = 0 ##TOTAL OF COLLIDED PACKETS
        nrProcessed = 0 ##TOTAL OF PROCESSED PACKETS
        nrReceived = 0 ###TOTAL OF RECEIVED PACKETS

if not mode_debbug:
    sys.stdout = old_stdout
    print("done LTb_"+str(p_skip_param))