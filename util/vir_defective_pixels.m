function [defective_pixels] = vir_defective_pixels()
defective_pixels = false(1,256,432);
defective_pixels(1, 15, 167      ) = true;
defective_pixels(1, 51, [212 280]) = true;
defective_pixels(1, 52, [280])     = true;
defective_pixels(1, 56, [174])     = true;
% defective_pixels(1, 66, [229:232]) = true;
defective_pixels(1, 81, [189:191]) = true;
defective_pixels(1, 82, [189:191]) = true;
defective_pixels(1, 83, [189:190]) = true;
defective_pixels(1, 85, 191)       = true;
defective_pixels(1, 87, 257)       = true;
defective_pixels(1, 96, 186:189)   = true;
defective_pixels(1, 97, 186:189)   = true;
defective_pixels(1, 98, 186:189)   = true;
defective_pixels(1,101, 223:225)   = true;
defective_pixels(1,102, 223:225)   = true;
defective_pixels(1,103, 223:225)   = true;
defective_pixels(1,120, 192:194)   = true;
defective_pixels(1,121, 192:194)   = true;
defective_pixels(1,122, 192:194)   = true;
defective_pixels(1,132, 182)       = true;
defective_pixels(1,139, 190:192)   = true;
defective_pixels(1,140, 191:192)   = true;
defective_pixels(1,149, 169:171)   = true;
defective_pixels(1,150, 169:171)   = true;
defective_pixels(1,169, 184)       = true;
defective_pixels(1,193, [245 246]) = true;
defective_pixels(1,198, [185 186]) = true;
defective_pixels(1,236, [186 187]) = true;
defective_pixels(1,245, [191 192]) = true;
defective_pixels(1,250, 291)       = true;

end