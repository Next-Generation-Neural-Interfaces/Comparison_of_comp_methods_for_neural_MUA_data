module clock_divier(
input clk,
input rst,
output reg clk_pos_2,
output reg clk_neg_2
);

always@(posedge clk or negedge rst) begin
	if(~rst)
		clk_pos_2 <= 0;
	else
		clk_pos_2 <= ~clk_pos_2;
end

always@(negedge clk or negedge rst) begin
	if(~rst)
		clk_neg_2 <= 0;
	else
		clk_neg_2 <= ~clk_neg_2;
end
endmodule


