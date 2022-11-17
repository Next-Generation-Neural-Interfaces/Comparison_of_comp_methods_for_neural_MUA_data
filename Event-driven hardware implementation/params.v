`define BIN_PERIOD 7// 50ms @ 7kHz
//`define BIN_PERIOD 35
//`define BIN_PERIOD 70
//`define BIN_PERIOD 140
//`define BIN_PERIOD 350
//`define BIN_PERIOD 700
`define BIN_PERIOD_WIDTH $clog2(`BIN_PERIOD-1)+1
`define SPIKE_NUMBER_WIDTH 2


`define MAX_CODEWORD_LENGTH 5
`define LENGTH_WIDTH 3


`define FREQ_BIT 6 //need to know that 
`define HISTOSIZE 64
`define HISTOCOUNTBIT 6 

`define ENCODER_LENGTH 4
`define TOTAL_LEN_BIT 16
`define ENCODER_NUM_BIT 2
`define ENOCDER_NUMBER 3



`define CH_NUM 1000
`define CH_BIT $clog2(`CH_NUM-1)+1


`define delta_channel 63
`define DCH_BIT $clog2(`delta_channel-1)+1

`define SPIKE_RATE_CLIP 5
`define SPIKE_RATE_BIT 3


`define buf_len `CH_NUM / `delta_channel
`define data_width 21