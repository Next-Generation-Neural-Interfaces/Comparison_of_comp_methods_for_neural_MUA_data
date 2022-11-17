`include "params.v"
module tb_delta_channel;
reg clk_,RST,cal;
wire clk, CLK;
reg detection[15000:0];
reg[14:0] counter,counter_;
wire detect;
wire[`SPIKE_RATE_BIT-1:0] spike_number,spike_number_in,spike_number_out,spike_number_out_,max_rate,max_rate_,rate_out;
wire [`FREQ_BIT-1:0] freq_count1,freq_count2,freq_count3,freq_count4,freq_count5,freq_1,freq_2,freq_3,freq_4,freq_5;
wire bin_finish;
//wire[`BIN_PERIOD_WIDTH-1:0] bin_period;
wire [`MAX_CODEWORD_LENGTH-1:0] codeword;
wire [`LENGTH_WIDTH-1:0] length;
wire finish_out ;

//wire [`SPIKE_RATE_BIT-1:0]spike_rate1,spike_rate2,spike_rate3,spike_rate4,spike_rate5;//spike_rate6,spike_rate7,spike_rate8,spike_rate9,spike_rate10,spike_rate11,spike_rate12,spike_rate13,spike_rate14,spike_rate15;
//wire [`FREQ_BIT-1:0] freq_count1,freq_count2,freq_count3,freq_count4,freq_count5,freq_count6,freq_count7,freq_count8,freq_count9,freq_count10,freq_count11,freq_count12,freq_count13,freq_count14,freq_count15;

wire [`SPIKE_RATE_BIT*2-1:0] ram_in, ram_out;
wire [`ENCODER_NUM_BIT-1:0] encoder_sel,encoder_sel_;

//wire [`SPIKE_RATE_BIT-1:0] max_rate;
wire [`CH_BIT-1:0] addr,channel_count;
reg [`CH_BIT-1:0] c_channel;
reg we;

initial begin
	$readmemb("binned_MUA_1_aligned_.txt",detection);
	clk_ <= 0;
	forever #5 clk_ <= ~clk_;
	
end

wire [0:`ENCODER_LENGTH*`SPIKE_RATE_CLIP-1] len;
always@(negedge CLK or negedge RST) //input counter after calibration
	if(~RST) counter <= 0;
	else 
		if (counter == 15000)
			counter <= 0;
		else 
			counter <= counter + 1; 

reg calibration,cali_finish;

wire triger;
assign detect = cali_finish?detection[counter_]:detection[counter];
assign triger = !cali_finish? finish_out:~CLK;
assign addr = !cali_finish? c_channel-`CH_BIT'b1:c_channel;
always@(negedge CLK or negedge RST) // input counter in calibration
	if(~RST||!cali_finish) counter_ <= 0;
	else 
		if (counter_ == 15000)
			counter_ <= 0;
		else 
			counter_ <= counter_ + 1; 
always@(posedge triger or negedge RST) begin //Channel counter after calibration
	if(~RST) begin
		c_channel <= 0;
	end
	else if(cali_finish == 0) begin
		c_channel <= c_channel + 1;
	end else 
		if(c_channel == `CH_NUM - 1 || c_channel == `CH_NUM )
			c_channel <= 0;
		else c_channel <= c_channel + 1;
end

always@(negedge finish_out or posedge calibration) begin //calibration signal logic. External calibration signal can cause the calibration and it finishes when all channels are calibrated.
	if(calibration)	cali_finish <= 0;
	else if(!finish_out && c_channel == `CH_NUM ) cali_finish <= 1;
	else cali_finish <= cali_finish;
end
	 	

 
binner_f BIN(CLK,RST,detect,cali_finish,spike_number_in,c_channel,spike_number,spike_number_out,spike_number_out_, bin_finish);

hist_unsort HIST( CLK&!cali_finish,RST,bin_finish,spike_number_out_,
freq_count1,freq_count2,freq_count3,freq_count4,freq_count5,max_rate, finish_out
);


assign ram_in = cali_finish?{spike_number, max_rate_}:{`SPIKE_RATE_BIT'b0, max_rate};
assign spike_number_in = ram_out[`SPIKE_RATE_BIT*2-1:`SPIKE_RATE_BIT];
assign max_rate_ = 3;
//assign encoder_sel_ = ram_out[`ENCODER_NUM_BIT-1:0];

dualRam DR(ram_in, cali_finish | finish_out,addr, ~clk, addr, clk, ram_out);

clock_divier CD(clk_,RST,clk,CLK);

reg[`SPIKE_RATE_BIT-1:0] temp_rate;
reg[`CH_NUM-1:0] channel_encoding = {`CH_NUM-2'b0,1'b1};
reg[`CH_BIT-1:0] delta_channel;
reg[`CH_BIT-1:0] last_channel;
wire not_max;
assign not_max = spike_number_out != max_rate_;
always@(negedge clk or negedge RST) begin //Event-driven
	if(!RST)	temp_rate <= 0;
	else if(not_max&&bin_finish)
		temp_rate <= spike_number_out;
	else temp_rate<=temp_rate;
end
always@(negedge clk or negedge RST) begin // delta channel calculation
	if(!RST) begin
		delta_channel <= 0;
		last_channel <= 0;
	end else if (!bin_finish)
		last_channel <= `CH_NUM-1;
	else if (not_max&&bin_finish) begin
		if(c_channel>last_channel) begin
			delta_channel <= c_channel - last_channel;
		end else
			delta_channel <= c_channel + `CH_NUM - last_channel; 
		last_channel <= c_channel;
	end
end

mapper MP(temp_rate,max_rate_,rate_out);


encoder ED(bin_finish,rate_out,codeword,length);

integer i,fd;
initial begin

//$monitor ("channel_count = %d READ = %b READY = %b DATA = %b", channel_count, mem_read, data_ready, data);
RST = 1;calibration = 0;
#1 RST = 0;calibration = 1; 
#1 RST = 1; calibration = 0;
//fd = $fopen("log.txt","w");
@(posedge cali_finish) $readmemb("binned_MUA_1_aligned.txt",detection);
fd = $fopen("log.txt","w");

//for(i = 0;i < 960 ;i = i+1) begin
@(posedge bin_finish)
@(posedge bin_finish)
for(i = 0;i < 96 ;i = i+1) begin
	@(negedge CLK)
		$fdisplay(fd, "time = %g channel_count = %d spike_number_out = %d max_rate = %d rate_out = %d codeword = %b length = %d \n",$time,channel_count,spike_number_out,max_rate,rate_out,codeword,length);
end

//@(negedge finish)  $display("%d %d",codeword,length);
//$display("%d:%b:%d",spike_number,codeword,length);
//@(negedge finish_out)  $display("%d",encoder_sel);
end
endmodule
  