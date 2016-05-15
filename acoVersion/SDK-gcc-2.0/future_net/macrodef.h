#ifndef MACRODEF_H
#define MACRODEF_H

#define BITX        0x0000
#define BIT0        0x0001
#define BIT1        0x0002
#define BIT2        0x0004
#define BIT3        0x0008
#define BIT4        0x0010
#define BIT5        0x0020
#define BIT6        0x0040
#define BIT7        0x0080
#define NODE_MAX_NUM          2000

#define NODE_MAX_ROUTS         100
#define INVALID                -1
#define NODE_NUMBER          BITX
#define NODE_TYPE            BIT0
#define SUB_NODE_INDEX      BITX
#define SUB_NODE_NEXT       BIT0
#define NODE_TYPE_COMM      BITX
#define NODE_TYPE_START     BIT0
#define NODE_TYPE_END       BIT1
#define NODE_TYPE_MID       BIT2


#define PHEROMONE_INIT            0.5
#define MAX_WEIGHT           	101
#define POSITIVE_CONTS            0.95
#define EVAPORATION_RATE          0.5
#define MIN_PHEROMONE             0.2
#define MAX_PHEROMONE             0.8
#define ANT_NUMBER                500
#define ANT_ALPHA                 1
#define ANT_BETA                  100
#define ANT_UNGO                  0
#define ANT_GOST                  1
#define ANT_TABU                  2
#define MAX_DIST_VAL              INT_MAX


#define hash_func(from,to)  (10000*from+to)
#endif
