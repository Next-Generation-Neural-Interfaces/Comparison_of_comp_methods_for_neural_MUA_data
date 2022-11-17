`include "params.v" 
/*Map the binned spike rate to the correct number for encodering */
module mapper(
input [`SPIKE_RATE_BIT-1:0] rate_in,
input [`SPIKE_RATE_BIT-1:0]max_rate,

output [`SPIKE_RATE_BIT-1:0] rate_out
);

wire [`SPIKE_RATE_BIT-1:0]spike_rate1,spike_rate2,spike_rate3,spike_rate4,spike_rate5;

assign spike_rate1 = max_rate == 0? 0:max_rate == 1? 1:max_rate == 2?3:max_rate == 3? 4:4;
assign spike_rate2 = max_rate == 0? 1:max_rate == 1? 0:max_rate == 2?1:max_rate == 3? 3:3;
assign spike_rate3 = max_rate == 0? 2:max_rate == 1? 2:max_rate == 2?0:max_rate == 3? 1:2;
assign spike_rate4 = max_rate == 0? 3:max_rate == 1? 3:max_rate == 2?2:max_rate == 3? 0:1;
assign spike_rate5 = max_rate == 0? 4:max_rate == 1? 4:max_rate == 2?4:max_rate == 3? 2:0;
//assign rate_out = rate_in == spike_rate1? 0:spike_rate2? 1:2;
assign rate_out = rate_in == spike_rate1? 0:rate_in ==spike_rate2? 1:rate_in ==spike_rate3?2:rate_in ==spike_rate4?3:4;
//assign rate_out = rate_in == spike_rate1? 0:spike_rate2? 1:spike_rate3?2:spike_rate4?3:spike_rate5?4:spike_rate6?5:6;
//assign rate_out = rate_in == spike_rate1? 0:spike_rate2? 1:spike_rate3?2:spike_rate4?3:spike_rate5?4:spike_rate6?5:spike_rate7?6:spike_rate8?7:8;
endmodule