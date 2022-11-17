`include "params.v"
/*
Binner that bins count the number of detections in certain period. It provides three different spike number signal for different purpose, which is comment in Line13-15
*/
module binner_f(
	input CLK,
	input RST,
	input detected,
	input cali_finish,
	input [`SPIKE_RATE_BIT-1:0] spike_number_in,
	//input [`BIN_PERIOD_WIDTH-1:0] bin_period,
	input  [`CH_BIT - 1:0] channel_count,
	output [`SPIKE_RATE_BIT-1:0]spike_number, // the spike number to be stored into RAM
	output reg [`SPIKE_RATE_BIT-1:0]spike_number_out, // the spike number to be compressed
	output reg [`SPIKE_RATE_BIT-1:0]spike_number_out_, // the spike number only used in calibration
	output finish
);

reg [`BIN_PERIOD_WIDTH-1:0] bin_count = 0;
assign finish = bin_count == `BIN_PERIOD-1;
assign spike_number =  finish? 0: spike_number_in == `SPIKE_RATE_CLIP-1? `SPIKE_RATE_CLIP-1 :  detected? spike_number_in + 1 : spike_number_in;

always@(posedge finish or negedge RST) begin
	if(~RST) begin
		spike_number_out_ <= 0;
	end else if(finish) begin
		if (detected) // need to be modified if spike number could overflow
			spike_number_out_ <= spike_number_out + 1;
		else
			spike_number_out_ <= spike_number_out;
	end else 
		spike_number_out_ <= spike_number_out_;
	
end

always@(posedge CLK or negedge RST) begin
	if(~RST) 
		spike_number_out <= 0;
	else if (!cali_finish) // during calibration, the spike number is count one channel after last channel finishes, so concurrently added spike number out
		if(finish)
			spike_number_out <= 0;
		else if (detected)
			spike_number_out <= spike_number_out + 1;
		else
			spike_number_out <= spike_number_out;
	else if(finish) begin // After calibration, the spike number is count channel by channel interchangably. So it adds 1 to the number reads from RAM (spike_number_in)
		if (detected) // need to be modified if spike number could overflow
			spike_number_out <= spike_number_in + 1;
		else
			spike_number_out <= spike_number_in;
	end else 
		spike_number_out <= spike_number_out;
	
end
always@(negedge CLK or negedge RST) begin
	if(~RST)
		bin_count <= 0;
	else if(!cali_finish) //during calibration, contiueously increase the count
		if(bin_count == `BIN_PERIOD-1)
			bin_count <= 0;
		else
			bin_count <= bin_count + 1;
	else if (channel_count == `CH_NUM - 1) // after calibration, only after all channels are counted, the bin_count is increased.
		if(bin_count == `BIN_PERIOD-1)
			bin_count <= 0;
		else
			bin_count <= bin_count + 1;
	else
		bin_count <= bin_count;
end

endmodule
