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
#import matplotlib.pyplot as plt
#import os
#import operator
#from mpl_toolkits import mplot3d
#from mpl_toolkits.mplot3d import Axes3D
#import PIL
import random
import re
import os
import datetime
import sys

name = "EB2"
mode_debbug = 0

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
    maxBSReceives = 500
    multi_nodes = [1500]
else:
    RANDOM_SEED = int(sys.argv[1])
    chan = int(sys.argv[2])
    packetlen = int(sys.argv[3])   ##NODES SEND PACKETS OF JUST 20 Bytes
    total_data = int(sys.argv[4]) ##TOTAL DATA ON BUFFER, FOR EACH NODE (IT'S THE BUFFER O DATA BEFORE START SENDING)
    beacon_time = int(sys.argv[5]) ###SAT SENDS BEACON EVERY CERTAIN TIME
    maxBSReceives = int(sys.argv[6]) ##MAX NUMBER OF PACKETS THAT BS (ie SATELLITE) CAN RECEIVE AT SAME TIME
    
    #multi_nodes = [int(sys.argv[7]), int(sys.argv[8]) ,int(sys.argv[9])]
    multi_nodes = [int(sys.argv[7]), int(sys.argv[8]) ,int(sys.argv[9]), int(sys.argv[10]),int(sys.argv[11]),int(sys.argv[12]),int(sys.argv[13]),int(sys.argv[14]),int(sys.argv[15]),int(sys.argv[16]),int(sys.argv[17]),int(sys.argv[18]),int(sys.argv[19]),int(sys.argv[20])]
random.seed(RANDOM_SEED) #RANDOM SEED IS FOR GENERATE ALWAYS THE SAME RANDOM NUMBERS (ie SAME RESULTS OF SIMULATION)
nodesToSend = []
packetsToSend = math.ceil(total_data/packetlen)
###GLOBAL PARAMS ####
bsId = 1 ##ID OF BASE STATION (NOT USED)

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
nrNoProcessed = 0 ##TOTAL OF INTRA-PACKETS NO PROCESSED
nrIntraTot = 0
nrLostMaxRec = 0
nrCollFullPacket = 0
nrSentIntra = 0 ##TOTAL OF SENT INTRA-PACKTES
nrReceivedIntra = 0 ##TOTAL OF RECEIVED INTRA-PACKETS


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

#DR = ["dr8","dr9","dr10","dr11"]
DR = [-137,-134.5,-134,-131.5]

## READ PARAMS FROM DIRECTORY ##
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

distance = np.concatenate((distance,distance,distance,distance,distance,distance,distance,distance,distance,distance,\
                           distance,distance,distance,distance,distance,distance,distance,distance,distance,distance,\
                           distance,distance,distance,distance,distance,distance,distance,distance,distance,distance,\
                           distance,distance,distance,distance,distance,distance,distance,distance,distance,distance,\
                           distance,distance,distance,distance,distance,distance,distance,distance,distance,distance,\
                           distance,distance,distance,distance,distance,distance,distance,distance,distance,distance,\
                           distance,distance,distance,distance,distance,distance,distance,distance,distance,distance,\
                           distance,distance,distance,distance,distance,distance,distance,distance,distance,distance,\
                           distance,distance,distance,distance,distance,distance,distance,distance,distance,distance,\
                           distance,distance,distance,distance,distance,distance,distance,distance,distance,distance,\
                           distance,distance,distance,distance,distance,distance,distance,distance,distance,distance,\
                           distance,distance,distance,distance,distance,distance,distance,distance,distance,distance,\
                           distance,distance,distance,distance,distance,distance,distance,distance,distance,distance,\
                           distance,distance,distance,distance,distance,distance,distance,distance,distance,distance,\
                           distance,distance,distance,distance,distance,distance,distance,distance,distance,distance))

Lpl = np.concatenate((Lpl,Lpl,Lpl,Lpl,Lpl,Lpl,Lpl,Lpl,Lpl,Lpl,\
                      Lpl,Lpl,Lpl,Lpl,Lpl,Lpl,Lpl,Lpl,Lpl,Lpl,\
                      Lpl,Lpl,Lpl,Lpl,Lpl,Lpl,Lpl,Lpl,Lpl,Lpl,\
                      Lpl,Lpl,Lpl,Lpl,Lpl,Lpl,Lpl,Lpl,Lpl,Lpl,\
                      Lpl,Lpl,Lpl,Lpl,Lpl,Lpl,Lpl,Lpl,Lpl,Lpl,\
                      Lpl,Lpl,Lpl,Lpl,Lpl,Lpl,Lpl,Lpl,Lpl,Lpl,\
                      Lpl,Lpl,Lpl,Lpl,Lpl,Lpl,Lpl,Lpl,Lpl,Lpl,\
                      Lpl,Lpl,Lpl,Lpl,Lpl,Lpl,Lpl,Lpl,Lpl,Lpl,\
                      Lpl,Lpl,Lpl,Lpl,Lpl,Lpl,Lpl,Lpl,Lpl,Lpl,\
                      Lpl,Lpl,Lpl,Lpl,Lpl,Lpl,Lpl,Lpl,Lpl,Lpl,\
                      Lpl,Lpl,Lpl,Lpl,Lpl,Lpl,Lpl,Lpl,Lpl,Lpl,\
                      Lpl,Lpl,Lpl,Lpl,Lpl,Lpl,Lpl,Lpl,Lpl,Lpl,\
                      Lpl,Lpl,Lpl,Lpl,Lpl,Lpl,Lpl,Lpl,Lpl,Lpl,\
                      Lpl,Lpl,Lpl,Lpl,Lpl,Lpl,Lpl,Lpl,Lpl,Lpl,\
                      Lpl,Lpl,Lpl,Lpl,Lpl,Lpl,Lpl,Lpl,Lpl,Lpl))

Prx = np.concatenate((Prx,Prx,Prx,Prx,Prx,Prx,Prx,Prx,Prx,Prx,\
                      Prx,Prx,Prx,Prx,Prx,Prx,Prx,Prx,Prx,Prx,\
                      Prx,Prx,Prx,Prx,Prx,Prx,Prx,Prx,Prx,Prx,\
                      Prx,Prx,Prx,Prx,Prx,Prx,Prx,Prx,Prx,Prx,\
                      Prx,Prx,Prx,Prx,Prx,Prx,Prx,Prx,Prx,Prx,\
                      Prx,Prx,Prx,Prx,Prx,Prx,Prx,Prx,Prx,Prx,\
                      Prx,Prx,Prx,Prx,Prx,Prx,Prx,Prx,Prx,Prx,\
                      Prx,Prx,Prx,Prx,Prx,Prx,Prx,Prx,Prx,Prx,\
                      Prx,Prx,Prx,Prx,Prx,Prx,Prx,Prx,Prx,Prx,\
                      Prx,Prx,Prx,Prx,Prx,Prx,Prx,Prx,Prx,Prx,\
                      Prx,Prx,Prx,Prx,Prx,Prx,Prx,Prx,Prx,Prx,\
                      Prx,Prx,Prx,Prx,Prx,Prx,Prx,Prx,Prx,Prx,\
                      Prx,Prx,Prx,Prx,Prx,Prx,Prx,Prx,Prx,Prx,\
                      Prx,Prx,Prx,Prx,Prx,Prx,Prx,Prx,Prx,Prx,\
                      Prx,Prx,Prx,Prx,Prx,Prx,Prx,Prx,Prx,Prx))

elev = np.degrees(np.arcsin(599/distance))

# =============================================================================
# IS7 = np.array([1,-8,-9,-9,-9,-9])
# IS8 = np.array([-11,1,-11,-12,-13,-13])
# IS9 = np.array([-15,-13,1,-13,-14,-15])
# IS10 = np.array([-19,-18,-17,1,-17,-18])
# IS11 = np.array([-22,-22,-21,-20,1,-20])
# IS12 = np.array([-25,-25,-25,-24,-23,1])
# IsoThresholds = np.array([IS7,IS8,IS9,IS10,IS11,IS12])
# 
# =============================================================================
ISDR8 = np.array([1,-23,-24,-25])
ISDR9 = np.array([-20,1,-20,-21])
ISDR10 = np.array([-18,-17,1,-17])
ISDR11 = np.array([-15,-14,-13,1])

IsoThresholds = np.array([ISDR8,ISDR9,ISDR10,ISDR11])

#THIS IS THE MATRIX OF CROSS INTERFERENCE BETWEEN SF

Collmap = [[0 for i in range(0,6)] for j in range(0,6)]

def simulate_scenario (nrNodes):
    env = simpy.Environment()
    
    def powerCollision_2(p1, p2):
        #powerThreshold = 6
        global Collmap
        #print ("SF: node {} has {} ; node {} has {}".format(p1.nodeid,p1.sf, p2.nodeid, p2.sf))
        #print ("pwr: node {} with rssi {} dBm and node {} with rssi {} dBm; diff {:3.2f} dBm".format(p1,p1.rssi[math.ceil(env.now)], p2.nodeid,p2.rssi[math.ceil(env.now)], p1.rssi[math.ceil(env.now)] - p2.rssi[math.ceil(env.now)]))
        if True: #p1.sf == p2.sf:
           if abs(p1.rssi[math.ceil(env.now)] - p2.rssi[math.ceil(env.now)]) < IsoThresholds[p1.dr-8][p2.dr-8]:
               print ("collision pwr both node {} and node {}".format(p1.nodeid, p2.nodeid))
               #Collmap[p1.sf-7][p2.sf-7] += 1
               #Collmap[p2.sf-7][p1.sf-7] += 1
               # packets are too close to each other, both collide
               # return both packets as casualties
               return (p1, p2)
           elif p1.rssi[math.ceil(env.now)] - p2.rssi[math.ceil(env.now)] < IsoThresholds[p1.dr-8][p2.dr-8]:
               # p2 overpowered p1, return p1 as casualty
               print ("collision pwr node {} overpowered node {}".format(p2.nodeid, p1.nodeid))
               print ("capture - p2 wins, p1 lost")
               #Collmap[p1.sf-7][p2.sf-7] += 1
               return (p1,)
           print ("capture - p1 wins, p2 lost")
           # p2 was the weaker packet, return it as a casualty
           #Collmap[p2.sf-7][p1.sf-7] += 1
           return (p2,)
               
    def checkcollision(packet):
        col = 0 # flag needed since there might be several collisions for packet
        processing = 0
        #print ("MAX RECEIVE IS: ", maxBSReceives)
        for i in range(0,len(packetsAtBS)):
            if packetsAtBS[i].header.processed == 1 or packetsAtBS[i].intraPacket.processed == 1:
                processing = processing + 1
        if (processing > maxBSReceives):
            packet.processed = 0
        else:
            packet.processed = 1
        if packetsAtBS:
            #print ("{:3.5f} || >> FOUND overlap... node {} (sf:{} bw:{} freq:{}) others: {}".format(env.now,packet.nodeid, packet.sf, packet.bw,packet.freq,len(packetsAtBS)))
            for other in packetsAtBS:
                if other.nodeid != packet.nodeid:
                   #print ("{:3.5f} || >> node {} overlapped with node {} (sf:{} bw:{} freq:{}). Let's check Freq...".format(env.now,packet.nodeid, other.nodeid, other.packet.sf, other.packet.bw,other.packet.freq))
                   # simple collision
                   #if frequencyCollision(packet, other.packet) and sfCollision(packet, other.packet):
                   if frequencyCollision(packet, other.header):# and timingCollision(packet, other.packet):
                       c = powerCollision_2(packet, other.header)
                       for p in c:
                          p.col = 1
                          if p == packet:
                             col = 1
                   if frequencyCollision(packet, other.intraPacket):# and timingCollision(packet, other.packet):
                       c = powerCollision_2(packet, other.intraPacket)
                       for p in c:
                          p.col = 1
                          if p == packet:
                             col = 1
            return col
        return 0
    
       
    def frequencyCollision(p1,p2):
        if (p1.ch == p2.ch):
            #print ("{:3.5f} || >> same channel for header on node {} and node {}.. Let's check sub-channels...".format(env.now,p1.nodeid,p2.nodeid))
            #if (p1.freqHopHeader[replica] == p2.freqHopHeader[replica]):
            if (p1.subCh == p2.subCh):
                #print ("{:3.5f} || >> same sub-channel for header on node {} and node {}".format(env.now,p1.nodeid,p2.nodeid))
                #print ("{:3.5f} || >> Header {} from node {} collided!!!".format(env.now,replica,p1.nodeid))
                return True
            else:
                #print ("{:3.5f} || >> No sub-channel collision".format(env.now))
                return False
        else:
            #print ("{:3.5f} || >> No header channel collision..".format(env.now))
            return False
    
    
    class myNode():
        def __init__(self, nodeid, bs, avgSendTime, packetlen, total_data):
            global channel
            global DR
            self.dr = 8
            #carriers = list(range(280))
            #random.shuffle(carriers) #TO CHOOSE THE HOPPING JUMPS
            self.nodeid = nodeid
            self.avgSendTime = avgSendTime
            self.bs = bs
            self.dist = distance[nodeid,:]
            self.elev = elev[nodeid,:]
            self.mindist = np.amin(distance[nodeid,:])
            self.mindist_pos = int(np.where(distance[nodeid,:] == np.amin(distance[nodeid,:]))[0])
            #print('node %d' %nodeid, "dist: ", self.dist[0])
            self.buffer = total_data
            self.packetlen = packetlen
            self.ch = int(random.choice(channel)) 
            self.packet = myPacket(self.nodeid, packetlen, self.dist)
            #self.freqHop = carriers[0:35]
            self.sent = 0 #INITIAL SENT PACKETS
            self.totalLost = 0 #INITIAL TOTAL LOST FOR PARTICULAR NODE
            self.totalColl = 0
            self.totalRec = 0
            self.totalProc = 0
            if self.dr == 8:
                carriers = list(range(280))
                random.shuffle(carriers) #TO CHOOSE THE HOPPING JUMPS
                self.freqHop = carriers[0:35]
            elif self.dr == 9:
                carriers = list(range(280))
                random.shuffle(carriers) #TO CHOOSE THE HOPPING JUMPS
                self.freqHop = carriers[0:35]
            elif self.dr == 10:
                carriers = list(range(688))
                random.shuffle(carriers) #TO CHOOSE THE HOPPING JUMPS
                self.freqHop = carriers[0:86]
            elif self.dr == 11:
                carriers = list(range(688))
                random.shuffle(carriers) #TO CHOOSE THE HOPPING JUMPS
                self.freqHop = carriers[0:86]
            
            self.header = myHeader(self.nodeid,self.dist,self.ch,self.freqHop, self.dr)
            self.intraPacket = myIntraPacket(self.nodeid,self.dist,self.ch,self.freqHop,self.dr)
            
            
    class myHeader ():
        def __init__(self,nodeid,dist,ch,freqHop,dr):
            global Ptx
            global Prx
            global Lpl
            global c
            global distance
            global channel
            global frequency
            self.nodeid = nodeid
            self.txpow = Ptx
            self.transRange = 150
            self.arriveTime = 0
            self.rssi = Prx[nodeid,:]
            self.rectime = 0.233
            #self.rectime = 1.5
            self.proptime = distance[nodeid,:]*(1/c)
            self.collided = 0
            self.noCollided = 0
            self.processed = 0
            self.noProcessed = 0
            self.ch = ch
            self.lost = bool
            self.Nlost = 0
            self.subCh = 0
            self.sentIntra = 0
            self.dr = dr
            self.col =0
            if dr == "dr8":
                self.freqHopHeader = freqHop[0:3]
            elif dr == "dr9":
                self.freqHopHeader = freqHop[0:2]
            elif dr == "dr10":
                self.freqHopHeader = freqHop[0:3]
            elif dr == "dr11":
                self.freqHopHeader = freqHop[0:2]
        
    
    class myIntraPacket ():
        def __init__(self,nodeid,dist,ch,freqHop,dr):
            global Ptx
            global Prx
            global Lpl
            global c
            global distance
            global channel
            global frequency
            self.nodeid = nodeid
            self.txpow = Ptx
            self.transRange = 150
            self.arriveTime = 0
            self.rssi = Prx[nodeid,:]
            self.freqHopIntraPacket = freqHop[3:]
            self.rectime = 50e-3
            #self.rectime = 3
            self.proptime = distance[nodeid,:]*(1/c)
            self.collided = 0
            self.noCollided = 0
            self.nrColl = 0
            self.processed = 0
            self.noProcessed = 0
            self.ch = ch
            self.lost = bool
            self.Nlost = 0
            self.subCh = 0
            self.sentIntra = 0
            self.dr = dr
            self.col =0
            if dr == "dr8":
                self.freqHopIntraPacket = freqHop[3:]
            elif dr == "dr9":
                self.freqHopIntraPacket = freqHop[2:]
            elif dr == "dr10":
                self.freqHopIntraPacket = freqHop[3:]
            elif dr == "dr11":
                self.freqHopIntraPacket = freqHop[2:]
    
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
            #print ("rectime node ", self.nodeid, "  ", self.rectime)
            #print ("Airtime for node {} is {} [seconds]".format(self.nodeid,self.rectime)) #from https://www.loratools.nl/#/airtime
            # denote if packet is collided
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
    
          
    def transmit(env,node):
        #while nodes[node.nodeid].buffer > 0.0:
        global wait_min
        global wait_max
        global back_off
        global beacon_time
        global logs
        global nodesToSend
        global DR
        while node.buffer > 0.0:
            #######STARTS TRANSMISSION AS DR8
            node.dr = 8
            node.sensi = -137
            carriers = list(range(280))
            random.shuffle(carriers) #TO CHOOSE THE HOPPING JUMPS
            node.freqHop = carriers[0:35]
            node.header.freqHopHeader = node.freqHop[0:3]
            node.intraPacket.freqHopIntraPacket = node.freqHop [3:]
            #######
            yield env.timeout(node.packet.rectime + float(node.packet.proptime[math.ceil(env.now)])) ##GIVE TIME TO RECEIVE BEACON
                          
            if node in packetsAtBS:
                print ("{:3.5f} || ERROR: packet is already in...".format(env.now))
            else:
                sensibility = sensi[12 - 7, [125,250,500].index(node.packet.bw) + 1]
                if node.packet.rssi[math.ceil(env.now)] < sensibility: #HERE WE ARE CONSIDERING RSSI AT TIME ENV.NOW
                    print ("{:3.5f} || Node {}: Can not reach beacon due Lpl".format(env.now,node.nodeid))
                    wait =0 ##LETS WAIT FOR NEXT BEACON
                    node.header.lost = False
                    node.intraPacket.lost = False
                    trySend = False
                    nIntraPackets = 0
                else:
                    nodesToSend.append(node.nodeid)
                    wait = random.uniform(1,back_off - node.packet.rectime - float(node.packet.proptime[math.ceil(env.now)])) ##TRIGGER BACK-OFF TIME
                    yield env.timeout(wait)
                    #print ("{:3.5f} || Node {} begins to transmit a packet".format(env.now,node.nodeid))
                    trySend = True
                    node.sent = node.sent + 1
                    node.buffer = node.buffer - node.packetlen
                    if node in packetsAtBS:
                        print ("{} || ERROR: packet is already in...".format(env.now))
                    else:
                        #sensibility = sensi[node.packet.sf - 7, [125,250,500].index(node.packet.bw) + 1]
                        sensibility = node.sensi
                        #print ("------Sensi is: ",sensibility)
                        if node.packet.rssi[math.ceil(env.now)] < sensibility: #HERE WE ARE CONSIDERING RSSI AT TIME ENV.NOW
                            print ("{:3.5f} || Node {}: The Packet will be Lost due Lpl".format(env.now,node.nodeid))
                            node.header.lost = True ## LOST ONLY CONSIDERING Lpl
                            node.intraPacket.lost = True ## LOST ONLY CONSIDERING Lpl
                            #nIntraPackets = 0
                            print ("###############lost !!!!!!!!")
                        else:
                            node.header.lost = False ## LOST ONLY CONSIDERING Lpl
                            node.intraPacket.lost = False ## LOST ONLY CONSIDERING Lpl
                            #print ("{:3.5f} || Prx for node {} is {:3.2f} dB".format(env.now, node.nodeid, node.packet.rssi[math.ceil(env.now)]))
                            #print ("Prx for node",node.nodeid, "is: ",node.packet.rssi[math.ceil(env.now)],"at time",env.now)
                           
                        for i in range(len(node.header.freqHopHeader)):
                            ###print ("{:3.5f} || Sending Header replica {} node {}...".format(env.now,i,node.nodeid))
                            ###print ("{:3.5f} || Let's try if there are collisions...".format(env.now))
                            node.header.subCh = node.header.freqHopHeader[i]
                            #print ("SUBCHANELLLL: ",node.header.subCh)
                            node.header.sentIntra +=1
                            isLost =0
                            if node.packet.rssi[math.ceil(env.now)] < sensibility:
                                node.header.Nlost +=1
                                isLost =1
                            if (checkcollision(node.header)==1):
                                #pass
                                if node.header.col == 1:
                                    if isLost == 0:
                                        node.header.collided +=1
                                #node.packet.collided = 1
                                #print ("---{:3.5f} || Collision for Header replica {} node {} !!!".format(env.now,i,node.nodeid))
                                #node.packet.collided = 1
                                #node.header.collided +=1 #ALREADY COUNTED IN FUNCTION                                
                            else:
                                ###print ("{:3.5f} || ...No Collision for Header replica {} node {}!".format(env.now,i,node.nodeid))
                                #node.packet.collided = 0
                                node.header.noCollided = 1 ##ALMOST ONE HEADER IS OK, THEN HEADER IS OK
                            packetsAtBS.append(node)
                            node.packet.addTime = env.now
                            isLost =0
                        
                            yield env.timeout(node.header.rectime)
                            if (node in packetsAtBS):
                                packetsAtBS.remove(node)
                        
                        ##CALCULATE N OF INTRAPACKETS BASED ON PACKETLEN
                        #payloadTime = airtime(12,1,node.packetlen,125)
                        if node.dr == 8 or node.dr == 10:
                            payloadTime = 1.85 - 0.233*3 
                        elif node.dr == 9 or node.dr ==11:
                            payloadTime = 1.07 - 0.233*2
                        
                        nIntraPackets = math.ceil(payloadTime / 50e-3)
                        #print ("NUMBER OF INTRA PACKETSSSS",nIntraPackets)
                        
                        for j in range (nIntraPackets):
                            ###print ("{:3.5f} || Sending intra-packet {} of {} for node {}...".format(env.now,j,nIntraPackets-1,node.nodeid))
                            ###print ("{:3.5f} || Let's try if there are collisions...".format(env.now))
                            node.intraPacket.subCh = node.intraPacket.freqHopIntraPacket[j]
                            node.intraPacket.sentIntra +=1
                            isLost =0
                            if node.packet.rssi[math.ceil(env.now)] < sensibility:
                                node.intraPacket.Nlost +=1
                                isLost =1
                            #print ("INTRA-PACKT SUB CHANNELLLL", node.intraPacket.subCh)
                            if (checkcollision(node.intraPacket)==1):
                                #pass
                                if node.intraPacket.col == 1:
                                    if isLost ==0:
                                        node.intraPacket.collided+=1
                                #print ("---{:3.5f} || Collision for intra-packet {} for node {} !!!".format(env.now,j,node.nodeid))
                                #node.intraPacket.collided+=1 #ALREADY COUNTED ON FUNCTION
                            else:
                                ###print ("{:3.5f} || ...No Collision for intra-packet {} for node {}!".format(env.now,j,node.nodeid))
                                node.intraPacket.noCollided +=1
                                pass
                            packetsAtBS.append(node)
                            node.packet.addTime = env.now
                            isLost =0
                            yield env.timeout(node.intraPacket.rectime)
                            if (node in packetsAtBS):
                                packetsAtBS.remove(node)
                            #print ("INTRA-PACKET NO-PROCESEDDD",node.intraPacket.noProcessed)
            
            node.header.noCollided = len(node.header.freqHopHeader)-node.header.Nlost-node.header.collided
            node.intraPacket.noCollided = nIntraPackets-node.intraPacket.Nlost-node.intraPacket.collided
            if node.header.noCollided <0:
                node.header.noCollided = 0
            if node.intraPacket.noCollided <0:
                node.intraPacket.noCollided = 0
                
            if trySend == 1:
                #print ("----count intra-packet collided", node.intraPacket.collided)
                if node.header.lost or node.intraPacket.lost:
                    logs.append("{:3.3f},{},{:3.3f},{:3.3f},{},PL,#{},#{},#{},#{}".format(env.now,node.nodeid,node.dist[math.ceil(env.now)],node.elev[math.ceil(env.now)],node.dr,nIntraPackets,node.intraPacket.noCollided,len(node.header.freqHopHeader),node.header.noCollided))
                else:
                    if node.dr ==8 or node.dr==10:
                        if node.header.collided == 3:
                            logs.append("{:3.3f},{},{:3.3f},{:3.3f},{},PCh,#{},#{},#{},#{}".format(env.now,node.nodeid,node.dist[math.ceil(env.now)],node.elev[math.ceil(env.now)],node.dr,nIntraPackets,node.intraPacket.noCollided,len(node.header.freqHopHeader),node.header.noCollided))
                        elif node.intraPacket.collided > (1/3)*nIntraPackets:
                            logs.append("{:3.3f},{},{:3.3f},{:3.3f},{},PCp,#{},#{},#{},#{}".format(env.now,node.nodeid,node.dist[math.ceil(env.now)],node.elev[math.ceil(env.now)],node.dr,nIntraPackets,node.intraPacket.noCollided,len(node.header.freqHopHeader),node.header.noCollided))
                        elif node.header.noProcessed == 3:
                            logs.append("{:3.3f},{},{:3.3f},{:3.3f},{},NP,#{},#{},#{},#{}".format(env.now,node.nodeid,node.dist[math.ceil(env.now)],node.elev[math.ceil(env.now)],node.dr,nIntraPackets,node.intraPacket.noCollided,len(node.header.freqHopHeader),node.header.noCollided))
                        elif node.intraPacket.noProcessed > (1/3)*nIntraPackets:
                            logs.append("{:3.3f},{},{:3.3f},{:3.3f},{},NP,#{},#{},#{},#{}".format(env.now,node.nodeid,node.dist[math.ceil(env.now)],node.elev[math.ceil(env.now)],node.dr,nIntraPackets,node.intraPacket.noCollided,len(node.header.freqHopHeader),node.header.noCollided))
                        else:
                            logs.append("{:3.3f},{},{:3.3f},{:3.3f},{},PE,#{},#{},#{},#{}".format(env.now,node.nodeid,node.dist[math.ceil(env.now)],node.elev[math.ceil(env.now)],node.dr,nIntraPackets,node.intraPacket.noCollided,len(node.header.freqHopHeader),node.header.noCollided))
                                   
                    elif node.dr==9 or node.dr==11:
                        if node.header.collided == 2:
                            logs.append("{:3.3f},{},{:3.3f},{:3.3f},{},PCh,#{},#{},#{},#{}".format(env.now,node.nodeid,node.dist[math.ceil(env.now)],node.elev[math.ceil(env.now)],node.dr,nIntraPackets,node.intraPacket.noCollided,len(node.header.freqHopHeader),node.header.noCollided))
                        elif node.intraPacket.collided > (2/3)*nIntraPackets:
                            logs.append("{:3.3f},{},{:3.3f},{:3.3f},{},PCp,#{},#{},#{},#{}".format(env.now,node.nodeid,node.dist[math.ceil(env.now)],node.elev[math.ceil(env.now)],node.dr,nIntraPackets,node.intraPacket.noCollided,len(node.header.freqHopHeader),node.header.noCollided))
                        elif node.header.noProcessed == 2:
                            logs.append("{:3.3f},{},{:3.3f},{:3.3f},{},NP,#{},#{},#{},#{}".format(env.now,node.nodeid,node.dist[math.ceil(env.now)],node.elev[math.ceil(env.now)],node.dr,nIntraPackets,node.intraPacket.noCollided,len(node.header.freqHopHeader),node.header.noCollided))
                        elif node.intraPacket.noProcessed > (2/3)*nIntraPackets:
                            logs.append("{:3.3f},{},{:3.3f},{:3.3f},{},NP,#{},#{},#{},#{}".format(env.now,node.nodeid,node.dist[math.ceil(env.now)],node.elev[math.ceil(env.now)],node.dr,nIntraPackets,node.intraPacket.noCollided,len(node.header.freqHopHeader),node.header.noCollided))
                        else:
                            logs.append("{:3.3f},{},{:3.3f},{:3.3f},{},PE,#{},#{},#{},#{}".format(env.now,node.nodeid,node.dist[math.ceil(env.now)],node.elev[math.ceil(env.now)],node.dr,nIntraPackets,node.intraPacket.noCollided,len(node.header.freqHopHeader),node.header.noCollided))
            
            ##RESET
            node.header.collided = 0
            node.header.processed = 0
            node.header.noProcessed = 0
            node.header.lost = False
            node.header.noCollided =0
            node.intraPacket.nrColl = 0
            node.intraPacket.collided = 0
            node.intraPacket.processed = 0
            node.intraPacket.noProcessed = 0
            node.intraPacket.lost = False
            node.intraPacket.noCollided = 0
            node.header.sentIntra = 0
            node.intraPacket.sentIntra = 0
            node.header.Nlost =0
            node.intraPacket.Nlost = 0
            
            if trySend:
                #print ("BEACON TIMEEE",beacon_time)
                #print ("WAITTT",wait)
                #print ("NODE HEADER TIME",node.header.rectime)
                #print ("ONE INTRA-PACKET TIMEE",node.intraPacket.rectime)
                #yield env.timeout(beacon_time-wait)
                yield env.timeout(beacon_time-wait-2*3*node.header.rectime-2*nIntraPackets*node.intraPacket.rectime)
            else:
                nIntraPackets = 0
                yield env.timeout(beacon_time-wait-3*node.header.rectime-nIntraPackets*node.intraPacket.rectime)
                    
    def beacon (env):
        global beacon_time
        global nodesToSend
        global logs
        i = 0
        while True:
            if i == 0:
                yield env.timeout(0)           
            else:
                yield env.timeout(beacon_time-2)
            i=i+1
            print ("{:3.5f} || ***A new beacon has been sended from Satellite***".format(env.now))
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
    return ([sent,nrCollFullPacket,None,None,nrReceived],logs)



#########################################################################
if chan == 1:
    ###SCENARIO 1 CHANNEL###
    channel = [0]
    
    nodes = [] ###EACH NODE WILL BE APPENDED TO THIS VARIABLE
    nrLost = 0 ### TOTAL OF LOST PACKETS DUE Lpl
    nrCollisions = 0 ##TOTAL OF COLLIDED PACKETS
    nrProcessed = 0 ##TOTAL OF PROCESSED PACKETS
    nrReceived = 0 ###TOTAL OF RECEIVED PACKETS
    nrNoProcessed = 0 ##TOTAL OF INTRA-PACKETS NO PROCESSED
    nrIntraTot = 0
    nrLostMaxRec = 0
    nrCollFullPacket = 0
    nrSentIntra = 0 ##TOTAL OF SENT INTRA-PACKTES
    nrReceivedIntra = 0 ##TOTAL OF RECEIVED INTRA-PACKETS
    
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
        nrNoProcessed = 0 ##TOTAL OF INTRA-PACKETS NO PROCESSED
        nrIntraTot = 0
        nrLostMaxRec = 0
        nrCollFullPacket = 0
        nrSentIntra = 0 ##TOTAL OF SENT INTRA-PACKTES
        nrReceivedIntra = 0 ##TOTAL OF RECEIVED INTRA-PACKETS


#########################################################################
if chan ==3:
    ###SCENARIO 3 CHANNELS###
    channel = [0,1,2]
    nodes = [] ###EACH NODE WILL BE APPENDED TO THIS VARIABLE
    nrLost = 0 ### TOTAL OF LOST PACKETS DUE Lpl
    nrCollisions = 0 ##TOTAL OF COLLIDED PACKETS
    nrProcessed = 0 ##TOTAL OF PROCESSED PACKETS
    nrReceived = 0 ###TOTAL OF RECEIVED PACKETS
    nrNoProcessed = 0 ##TOTAL OF INTRA-PACKETS NO PROCESSED
    nrIntraTot = 0
    nrLostMaxRec = 0
    nrCollFullPacket = 0
    nrSentIntra = 0 ##TOTAL OF SENT INTRA-PACKTES
    nrReceivedIntra = 0 ##TOTAL OF RECEIVED INTRA-PACKETS
    i =0
    scenario_3ch = np.zeros((len(multi_nodes),5))
    results = []
    
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
        nrNoProcessed = 0 ##TOTAL OF INTRA-PACKETS NO PROCESSED
        nrIntraTot = 0
        nrLostMaxRec = 0
        nrCollFullPacket = 0
        nrSentIntra = 0 ##TOTAL OF SENT INTRA-PACKTES
        nrReceivedIntra = 0 ##TOTAL OF RECEIVED INTRA-PACKETS

if not mode_debbug:
    sys.stdout = old_stdout
    print("done ",name)





