### Version v2.5


import sys
import os.path
import argparse
import pandas as pd
import numpy as np
from scipy import interpolate
from binarytree import Node
import svgwrite
import seaborn as sns


def df2dict(data_df,exp_col):
	data_dict = dict()

	exps = data_df[exp_col].groupby(data_df['cell'])

	for cell_name,inds in exps.groups.items():
		exp_list = list(data_df[exp_col][inds])
		start_time = data_df['time'][inds[0]]
		time_length = len(exp_list)
		exp_avg = sum(exp_list)/time_length

		data_dict[cell_name] = {'stp':start_time,'tp_num':time_length,'exp_avg':exp_avg,'exp_list':exp_list}

	return data_dict


def to_tree_dict(data_dict):
	data_tree_dict = dict()
	root = tree_init()
	for cell in root.levelorder:
		data_tree_dict[cell.name] = cell

	### Build node for every cell
	for cell_name, cell_data in data_dict.items():
		data_tree_dict[cell_name] = build_node(cell_name,cell_data)

	### Relate the parent and the child
	for cell_name in data_dict:
		cell_parent_name,sub_direction = find_parent(cell_name)

		if sub_direction:
			data_tree_dict[cell_parent_name].right = data_tree_dict[cell_name]
		else:
			data_tree_dict[cell_parent_name].left = data_tree_dict[cell_name]

	return data_tree_dict


def tree_init():

	init_tp = 5

	root = Node(init_tp)
	root.name = 'P0'
	root.exp_list = [0]*init_tp
	root.stp= 0

	root.left = Node(init_tp)
	root.left.name = 'AB'
	root.left.exp_list = [0]*init_tp
	root.left.stp= 0

	root.right = Node(init_tp)
	root.right.name = 'P1'
	root.right.exp_list = [0]*init_tp
	root.right.stp= 0

	root.left.left = Node(init_tp)
	root.left.left.name = 'ABa'
	root.left.left.exp_list = [0]*init_tp
	root.left.left.stp= 0

	root.left.right = Node(init_tp)
	root.left.right.name = 'ABp'
	root.left.right.exp_list = [0]*init_tp
	root.left.right.stp= 0

	root.right.left = Node(init_tp)
	root.right.left.name = 'EMS'
	root.right.left.exp_list = [0]*init_tp
	root.right.left.stp= 0

	root.right.right = Node(init_tp)
	root.right.right.name = 'P2'
	root.right.right.exp_list = [0]*init_tp
	root.right.right.stp= 0

	return root


def build_node(cell_name,cell_data):
	a_node = Node(cell_data['tp_num'])
	a_node.name = cell_name
	a_node.stp = cell_data['stp']
	a_node.exp_avg = cell_data['exp_avg']
	a_node.exp_list = cell_data['exp_list']
	return a_node


def find_parent(cell_name):
	if cell_name in ('Z2','Z3'):
		return 'P4',('Z2','Z3').index(cell_name)
	elif cell_name in ('D','P4'):
		return 'P3',('D','P4').index(cell_name)
	elif cell_name in ('C','P3'):
		return 'P2',('C','P3').index(cell_name)
	elif cell_name in ('MS','E'):
		return 'EMS',('MS','E').index(cell_name)
	elif cell_name in ('EMS','P2'):
		return 'P1',('EMS','P2').index(cell_name)
	elif cell_name in ('AB','P1'):
		return 'P0',('AB','P1').index(cell_name)
	elif cell_name == 'P0':
		return None,None
	else:
		if cell_name[-1] in ('a','p'):
			return cell_name[:-1],('a','p').index(cell_name[-1])
		elif cell_name[-1] in ('l','r'):
			return cell_name[:-1],('l','r').index(cell_name[-1])
		elif cell_name[-1] in ('d','v'):
			return cell_name[:-1],('d','v').index(cell_name[-1])


def get_endtp(tree_dict):

	endtp = 0

	for cell in tree_dict['P0'].levelorder:
		if cell.stp > 0:
			endtp_new = cell.stp + cell.value - 1
			if endtp < endtp_new:
				endtp = endtp_new

	return endtp


def get_tree_value_depth(tree_dict,root):

	value_depth = dict()

	for cell in tree_dict[root].levelorder:
		value_depth[cell.name] = len(cell.exp_list)
		if cell.name != root:
			value_depth[cell.name] += value_depth[find_parent(cell.name)[0]]

	return max(value_depth.values())


def interpolation(ref,align):
	for cell in ref['P0'].levelorder:
		if cell.name in align.keys():
			ref_exp_list = cell.exp_list
			align_exp_list = align[cell.name].exp_list
			ref_len = len(ref_exp_list)
			align_len = len(align_exp_list)
			if cell.left and cell.right:
				if ref_len < 2:
					align_exp_list = align_exp_list[:ref_len]
				elif align_len < 2:
					align_exp_list = align_exp_list + [align_exp_list[-1]] * (ref_len - align_len)
				elif ref_len != align_len:
					x = np.linspace(0,1,len(align_exp_list))
					y = align_exp_list
					f = interpolate.interp1d(x,y,kind='linear') ### linear, nearest, nearest-up, zero, slinear, quadratic, cubic, previous, next
					xnew = np.linspace(0,1,len(ref_exp_list))
					ynew = f(xnew)
					align_exp_list = ynew
			else:
				if ref_len < align_len:
					align_exp_list = align_exp_list[:ref_len]
				elif ref_len > align_len:
					align_exp_list = align_exp_list + [align_exp_list[-1]] * (ref_len - align_len)
			
			align[cell.name].exp_list = align_exp_list
	return align


def get_linecolors():
	clvl = args.clvl
	clvlmax = args.clvlmax
	clvllow = args.clvllow
	clvllowplus = args.clvllowplus
	clvlhigh = args.clvlhigh

	assert clvl > 4, 'The color level should be larger than 4.'

	if clvlmax:
		assert clvlmax > 4, 'The max color level should be larger than 4.'
		assert clvlmax >= clvl, 'The clvlmax should not be smaller than clvl.'
	else:
		clvlmax = clvl

	if clvlhigh == 0:
		clvlhigh = clvlmax-(clvl-clvllowplus)+1

	assert clvlhigh <= clvlmax-(clvl-clvllowplus)+1, 'The clvlhigh is invalid. The clvlhigh should not be larger than {}'.format(clvlmax-(clvl-clvllowplus)+1) 

	color_set_max = list()
	color_set = list()

	if args.seaborn:
		if args.revc:
			color_set_max.append(sns.dark_palette((0,args.cl1/255,0),clvlmax))
			color_set_max.append(sns.dark_palette((args.cl1/255,0,0),clvlmax))
		else:
			color_set_max.append(sns.dark_palette((args.cl1/255,0,0),clvlmax))
			color_set_max.append(sns.dark_palette((0,args.cl1/255,0),clvlmax))
	else:
		if args.revc:
			color_set_max.append([(0.0,i,0.0) for i in np.linspace(0,args.cl1/255,clvlmax)])
			color_set_max.append([(i,0.0,0.0) for i in np.linspace(0,args.cl2/255,clvlmax)])
		else:
			color_set_max.append([(i,0.0,0.0) for i in np.linspace(0,args.cl1/255,clvlmax)])
			color_set_max.append([(0.0,i,0.0) for i in np.linspace(0,args.cl2/255,clvlmax)])
	
	for i in range(len(color_set_max)):
		color_set.append(list())
		color_set[i].extend(color_set_max[i][(clvllow-1):(clvllow-1+clvllowplus)])
		color_set[i].extend(color_set_max[i][(clvlhigh-1):(clvlhigh-1+(clvl-clvllowplus))])

	return color_set

def set_linecolor(exp_i,color_palette):
	clvl = args.clvl
	low_exp = 0
	high_exp = 5000

	rate = (exp_i - low_exp)/(high_exp - low_exp)

	if rate >= 1:
		color = color_palette[-1]
	elif rate <= 0:
		color = color_palette[0]
	else:
		color = color_palette[int((exp_i - low_exp)/(high_exp - low_exp)*(clvl-2))+1]

	return(color)


def blend_color(color1, color2):
	r1, g1, b1 = color1
	r2, g2, b2 = color2

	r = max(r1,r2)
	g = max(g1,g2)
	b = max(b1,b2)

	return(r, g, b)


def parameters_setting(ref,root):
	(marginleft, marginright, margintop, marginbottom) = (100, 100, 100, 100)
	(axis_width,brand_width) = (80,80)

	linewidth = args.linewidth
	lineinter = args.lineinter
	scale = 5

	width = marginleft + marginright + (ref[root].leaf_count)*(linewidth+lineinter) - lineinter + axis_width + brand_width
	height = margintop + marginbottom + scale*get_tree_value_depth(ref,root)


	xleaf = marginleft + axis_width
	for cell in ref[root].postorder:
		if cell.left and cell.right:
			cell.posx = cell.left.posx + 0.5*(cell.right.posx - cell.left.posx)
		elif not (cell.left or cell.right):
			cell.posx = xleaf
			xleaf += lineinter + linewidth
		else:
			print('The data of cell {} is not valid!'.format(cell.name))

	label_color = svgwrite.rgb(*[160/255*100]*3, '%')

	return (marginleft, marginright, margintop, marginbottom, axis_width, brand_width, scale, width, height, linewidth, ref, label_color)




def draw_tree(root, filenames, ref, align=None):
	def draw_axis():
		endtp = get_endtp(ref)
		starttp = ref[root].stp

		axis_x1 = marginleft
		axis_x2 = axis_x1 + 10
		axis_y1 = margintop

		if ref[root].stp == 0:
			for cell in ref[root].levelorder:
				if cell.stp > 0:
					axis_y1 = cell.posy
					starttp = cell.stp
					break

		axis_y4 = height - marginbottom
		axis_y2 = (axis_y4 - axis_y1)/3 + axis_y1
		axis_y3 = 2*(axis_y4 - axis_y1)/3 + axis_y1
		axis_x0 = axis_x1 - 10
		axis_y5 = (axis_y4 - axis_y1)/2 + axis_y1 + 50
		dwg.add(dwg.line((axis_x1, axis_y1), (axis_x1, axis_y4), stroke='black', stroke_width="1px"))

		dwg.add(dwg.line((axis_x1, axis_y1), (axis_x2, axis_y1), stroke='black', stroke_width="1px"))
		dwg.add(dwg.line((axis_x1, axis_y2), (axis_x2, axis_y2), stroke='black', stroke_width="1px"))
		dwg.add(dwg.line((axis_x1, axis_y3), (axis_x2, axis_y3), stroke='black', stroke_width="1px"))
		dwg.add(dwg.line((axis_x1, axis_y4), (axis_x2, axis_y4), stroke='black', stroke_width="1px"))

		x_shift,y_shift = 5,7
		dwg.add(dwg.text(str(starttp), insert=(axis_x2+x_shift, axis_y1+y_shift), fill='black', style = "font-size:25px; font-family:Times New Roman"))
		dwg.add(dwg.text(str(int(starttp+(endtp-starttp)/3)), insert=(axis_x2+x_shift, axis_y2+y_shift), fill='black', style = "font-size:25px; font-family:Times New Roman"))
		dwg.add(dwg.text(str(int(starttp+(endtp-starttp)/3*2)), insert=(axis_x2+x_shift, axis_y3+y_shift), fill='black', style = "font-size:25px; font-family:Times New Roman"))
		dwg.add(dwg.text(str(endtp), insert=(axis_x2+x_shift, axis_y4+y_shift), fill='black', style = "font-size:25px; font-family:Times New Roman"))

		dwg.add(dwg.text('Timepoint', insert=(axis_x0, axis_y5), fill='black', transform = "rotate(-90,{},{})".format(axis_x0,axis_y5), style = "font-size:25px; font-family:Times New Roman"))

	def draw_brand():
		brand_x = width - marginright - brand_width + 50
		brand_ys = margintop
		brand_yn = height - marginbottom
		brand_y_step = (brand_yn - brand_ys)/len(color_set[0])
		brand_ym = (brand_yn - brand_ys)/2 + brand_ys - 90


		for i,c_set in enumerate(color_set):
			brand_x += i*15
			brand_y1 = brand_ys
			for color in c_set:
				color = svgwrite.rgb(*[i*100 for i in color], '%')

				brand_y2 = brand_y1 + brand_y_step

				dwg.add(dwg.line((brand_x, brand_y1), (brand_x, brand_y2), stroke=color, stroke_width="15px"))

				brand_y1 = brand_y2

		brand_x += 15
		dwg.add(dwg.text('Expression Level', insert=(brand_x, brand_ym), fill='black', transform = "rotate(90,{},{})".format(brand_x, brand_ym), style = "font-size:25px; font-family:Times New Roman"))

	def draw_title():
		dwg.add(dwg.text('_'.join(filenames), insert=(10, 30), fill='black', style = "font-size:25px; font-family:Times New Roman;font-style:italic"))
		dwg.add(dwg.text('{}-cell stage'.format(ref['P0'].leaf_count), insert=(10, 30 + 30), fill='black', style = "font-size:25px; font-family:Times New Roman;font-style:italic"))

	color_set = get_linecolors()

	if align:
		align = interpolation(ref,align)
	else:
		del color_set[-1]

	(marginleft, marginright, margintop, marginbottom, axis_width, brand_width, scale, width, height, linewidth, ref, label_color) = parameters_setting(ref,root)


	leader_cells = ['P0','P1','P2','P3','P4','AB','ABa','ABp','EMS','E','MS','C','D','Z2','Z3']

	dwg = svgwrite.Drawing(size = ("{}px".format(width), "{}px".format(height)))

	ref[root].posy = margintop

	for cell in ref[root].levelorder:
		x1 = cell.posx
		x2 = x1
		y1 = cell.posy
		y1_copy = y1
		y2 = y1


		for i,exp_i in enumerate(cell.exp_list):
			y2 = y2+1*scale
			color = set_linecolor(exp_i,color_set[0])

			if align and cell.name in align.keys():
				color2 = set_linecolor(align[cell.name].exp_list[i],color_set[1])
				color = blend_color(color, color2)

			color = svgwrite.rgb(*[i*100 for i in color], '%')
			dwg.add(dwg.line((x1, y1), (x2, y2), stroke=color, stroke_width="{}px".format(linewidth)))

			y1 = y2

		if cell.left and cell.right:
			dwg.add(dwg.line((cell.left.posx, y2), (cell.right.posx, y2), stroke='black', stroke_width="2px"))

			x_label = x1 + linewidth
			y_label = y1_copy + scale

			if args.label and cell.name not in leader_cells and cell.name != root:
				dwg.add(dwg.text(cell.name, insert=(x_label, y_label), fill=label_color, transform = "rotate(90,{},{})".format(x_label,y_label), style = "font-size:13px; font-family:Times New Roman"))

			cell.left.posy = y2
			cell.right.posy = y2
		elif args.label:
			x_label = x1 - linewidth
			y_label = height - marginbottom + 2*scale
			dwg.add(dwg.text(cell.name, insert=(x_label, y_label), fill='black', transform = "rotate(45,{},{})".format(x_label,y_label), style = "font-size:10px; font-family:Times New Roman"))


		if cell.name == root:
			y_label = y1_copy - scale
			x_label = x1 - 12.5
			dwg.add(dwg.text(cell.name, insert=(x_label, y_label), fill='black', style = "font-size:25px; font-family:Times New Roman"))
		elif cell.name in leader_cells:
			y_label = y1_copy
			cell_parent_name,sub_direction = find_parent(cell.name)
			if sub_direction:
				x_label = x1 + 5
				dwg.add(dwg.text(cell.name, insert=(x_label, y_label), fill='black', style = "font-size:18px; font-family:Times New Roman"))
			else:
				x_label = x1 - 13*len(cell.name)
				dwg.add(dwg.text(cell.name, insert=(x_label, y_label), fill='black', style = "font-size:18px; font-family:Times New Roman"))



	if args.axis:
		draw_axis()

	if args.brand:
		draw_brand()

	if args.title:
		draw_title()


	dwg.saveas("{}/".format(args.output_file)+'_'.join(filenames+[args.exp_col,str(get_endtp(ref))+'t',str(ref['P0'].leaf_count)+'c',root])+'.svg')





def main():


	endtp = args.endtp
	cellstage = args.cellstage

	filenames = [os.path.splitext(os.path.basename(args.ref_file))[0]]


	### Load data
	print("=====> Using '{}' column.".format(args.exp_col))
	series = list()

	### Load reference data
	print("=====> Loading reference data ...")
	print('  -> reference: {}'.format(args.ref_file))

	ref_df = pd.read_csv(args.ref_file)

	### Trim reference data to be draw
	if cellstage:
		leafnum = 0
		endtp = 130
		endtp_record = dict()
		while leafnum != cellstage:

			if endtp in endtp_record.keys():
				break

			ref_df_trim = ref_df[ref_df['time']<=endtp]

			### DataFrame to Dict
			ref_dict = df2dict(ref_df_trim,args.exp_col)

			### Dict to Tree Dict
			ref_tree_dict = to_tree_dict(ref_dict)

			leafnum = ref_tree_dict['P0'].leaf_count

			endtp_record[endtp]=abs(leafnum - cellstage)
			print('     timepoint:{}, cellstage:{}'.format(endtp,leafnum))

			if leafnum < cellstage:
				endtp += 1
			elif leafnum > cellstage:
				endtp -= 1

		endtp = min(endtp_record.keys(),key=lambda x: endtp_record[x])

	if endtp:
		ref_df_trim = ref_df[ref_df['time']<=endtp]
	else:
		ref_df_trim = ref_df

	### DataFrame to Dict
	ref_dict = df2dict(ref_df_trim,args.exp_col)
	### Dict to Tree Dict
	ref_tree_dict = to_tree_dict(ref_dict)
	series.append(ref_tree_dict)
	print('     draw(timepoint:{}, cellstage:{})'.format(get_endtp(ref_tree_dict),ref_tree_dict['P0'].leaf_count))

	### Load aligns data
	if args.align_file:
		filenames.append(os.path.splitext(os.path.basename(args.align_file))[0])
		print("=====> Loading align data ...")
		print('  -> aligns: {}'.format(args.align_file))
		align_df = pd.read_csv(args.align_file)
		### DataFrame to Dict
		align_dict = df2dict(align_df,args.exp_col)
		### Dict to Tree Dict
		align_tree_dict = to_tree_dict(align_dict)
		series.append(align_tree_dict)
		print('     draw(timepoint:{}, cellstage:{})'.format(get_endtp(align_tree_dict),align_tree_dict['P0'].leaf_count))


	print("=====> Drawing tree...")
	draw_tree(args.root, filenames, *series)


if __name__ == '__main__':

	parser = argparse.ArgumentParser(description='Draw cell lineage tree!')
	parser.add_argument("-c", "--column", type=str, dest="exp_col", default='blot', help="Select a specific column to show on the tree")
	parser.add_argument("-r", "--root", type=str, dest="root", default='P0', help="Specify the root of the tree")
	parser.add_argument("--ref", type=str, dest="ref_file", help="Reference CD file")
	parser.add_argument("--align", type=str, dest="align_file", help="Align CD file")
	parser.add_argument("-o", "--output", type=str, dest="output_file", default='results', help="Output file folder")
	parser.add_argument("-nl", "--nolabel", action='store_false', dest="label", default=1, help="Draw cell label")
	parser.add_argument("-na","--noaxis", action='store_false', dest="axis", default=1, help="Draw timepoint axis")
	parser.add_argument("-nb","--nobrand", action='store_false', dest="brand", default=1, help="Draw brand")
	parser.add_argument("-nt","--notitle", action='store_false', dest="title", default=1, help="Draw title")
	group = parser.add_mutually_exclusive_group()
	group.add_argument("--endtp", type=int, dest="endtp", help="Set an end time point of the tree")
	group.add_argument("--cellstage", type=int, dest="cellstage", help="Set a cellstage of the tree")
	parser.add_argument("--linewidth", type=int, dest="linewidth", default=5, help="Set the linewidth")
	parser.add_argument("--lineinter", type=int, dest="lineinter", default=15, help="Set the lineinter")
	parser.add_argument("--clvl", type=int, dest="clvl", default=10, help="The number of color levels to be draw")
	parser.add_argument("--clvlmax", type=int, dest="clvlmax", default=0, help="The number of color levels in the palette")
	parser.add_argument("--clvllow", type=int, dest="clvllow", default=1, help="Low brightness part start position")
	parser.add_argument("--clvllowplus", type=int, dest="clvllowplus", default=1, help="The number of levels in the low brightness part")
	parser.add_argument("--clvlhigh", type=int, dest="clvlhigh", default=0, help="High brightness part start position")
	parser.add_argument("--revc", action='store_true', dest="revc", default=0, help="Reverse the reference and the align color")
	parser.add_argument("--cl1", type=int, dest="cl1", default=255, help="The first color max brightness")
	parser.add_argument("--cl2", type=int, dest="cl2", default=255, help="The second color max brightness")
	parser.add_argument("--seaborn", action='store_true', dest="seaborn", default=0, help="Use seaborn color palette")

	args = parser.parse_args()


	main()


