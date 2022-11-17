`include "params.v"
module dualRam (din, write_en, waddr, wclk, raddr, rclk, dout);//512x8 
 parameter addr_width = `CH_BIT;
 parameter data_width = `SPIKE_RATE_BIT*2 ;
 input [addr_width-1:0] waddr, raddr; 
 input [data_width-1:0] din;
 input write_en, wclk, rclk; 
 output reg [data_width-1:0] dout;
 reg [data_width-1:0] mem [0:`CH_NUM-1]/* synthesis syn_ramstyle = "block_ram" */;
initial 
$readmemh ("mem.txt",mem);

always @(posedge wclk) // Write memory. 
 begin
 if (write_en)
 mem[waddr] <= din; // Using write address bus.
 end
 always @(posedge rclk) // Read memory. 
 begin
 dout <= mem[raddr]; // Using read address bus. 
 end
endmodule