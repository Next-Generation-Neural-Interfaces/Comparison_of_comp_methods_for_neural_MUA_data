`include "params.v"
/*A FSM count the histogram of binning result. Records the spike number with highest frequency*/
module hist_unsort(
input CLK,
input RST,
input finish_in,
input [`SPIKE_RATE_BIT-1:0] rate,


output [`FREQ_BIT-1:0] freq_count1,
output [`FREQ_BIT-1:0] freq_count2,
output [`FREQ_BIT-1:0] freq_count3,
output [`FREQ_BIT-1:0] freq_count4,
output [`FREQ_BIT-1:0] freq_count5,
output reg [`SPIKE_RATE_BIT-1:0] max_rate,

output reg finish
);
wire [`SPIKE_RATE_BIT-1:0] rate_in;
wire internal_rst;
reg rst = 1;


reg [`FREQ_BIT-1:0] freq_count[0:`SPIKE_RATE_CLIP];
wire [`FREQ_BIT-1:0] freq_temp = rate_in == 0? freq_count[0]:rate_in == 1? freq_count[1]:rate_in == 2? freq_count[2]:rate_in == 3? freq_count[3]:freq_count[4];  //  These two lines can reduced the number of logic for access these numbers.
wire [`FREQ_BIT-1:0] freq_max = max_rate == 0? freq_count[0]:max_rate == 1? freq_count[1]:max_rate == 2? freq_count[2]:max_rate == 3? freq_count[3]:freq_count[4];

reg [`HISTOCOUNTBIT-1:0] histo_count = 0;

assign internal_rst = rst & RST;

assign finish_out = histo_count == `HISTOSIZE-1;

assign rate_in = rate>`SPIKE_RATE_CLIP? `SPIKE_RATE_CLIP:rate;
parameter IDLE = 3'd0;
parameter PLUS1 = 3'd1;
parameter OVFLOW = 3'd2;
parameter RESET = 3'd3;



reg[2:0] st_next;
reg[2:0] st_cur;

always@(posedge CLK or negedge RST) begin
	if (~RST) begin
            st_cur      <= 'b0 ;
	end
        else
            st_cur      <= st_next ;
end

always@(*) begin
	st_next = st_cur;
	case(st_cur)
		IDLE:
			case(finish_in) //is binning finished?

			1:	st_next = PLUS1;

			0:	st_next = IDLE;
			   
			endcase
		PLUS1:	case(finish_out) // is histogram full?
			1: st_next = OVFLOW;
			0: st_next = IDLE;
			endcase
		OVFLOW: st_next = RESET;
		RESET: st_next = IDLE;
	endcase
end
			
always@(posedge finish_in or negedge RST) begin
	if(!RST||histo_count == `HISTOSIZE-1)
		histo_count <= 0;
	else
		histo_count <= histo_count + 1;
end
integer i;
always@(posedge CLK or negedge RST) begin
	if(!RST) begin	
		for (i = 0; i<=`SPIKE_RATE_CLIP; i = i + 1) begin
			freq_count[i] <= 0;
		end
		max_rate <= 0;
	end
	else begin
		if(st_cur == PLUS1)begin // if binning finish, increase the related spike rate frequency and record the currently highest frequency rate
			freq_count[rate_in] <= freq_temp+ 1;
			if (freq_max<freq_temp + 1)
				max_rate <= rate_in;
		end else if (st_cur == OVFLOW)	
			finish <= 1;
		else if(st_cur == RESET) begin
			finish <= 0;
			for (i = 0; i<=`SPIKE_RATE_CLIP; i = i + 1)
				freq_count[i] <= 0;			
			max_rate <= 0;
		end 
	end

			
end

assign freq_count1 = freq_count[0];
assign freq_count2 = freq_count[1];
assign freq_count3 = freq_count[2];
assign freq_count4 = freq_count[3];
assign freq_count5 = freq_count[4];


endmodule