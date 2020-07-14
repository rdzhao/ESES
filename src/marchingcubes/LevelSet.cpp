#include "LevelSet.h"

int sign(double num)
{
	if (num == fabs(num))
		return 1;
	else
		return -1;
}

int input_levelset_grid_info(std::string filename, MarchingCubes& mc)
{
	//read
	std::fstream file;
	file.open(filename.c_str(), std::ios::in);
	std::string line;
	for (int i = 0; i < 6; i++)
		getline(file, line);
	getline(file, line);
	std::stringstream ss(line);
	std::string trash;
	for(int i=0; i<5; i++)
	  ss>>trash;
	
	int new_size_x, new_size_y, new_size_z;
	ss>>new_size_x;
	ss>>new_size_y;
	ss>>new_size_z;
	getline(file, line);
	std::vector<double> grid;
	//getline(file, line);
	//std::stringstream ss(line);
	//string token_d;
	//ss>>token_d;
	//cout<<token_d<<endl;
	
	while (getline(file, line))
	{
		std::stringstream ss(line);
		double token_d;
		ss >> token_d;
		grid.push_back(token_d);
	}
	
	//cout<<"test: "<<grid[102432]<<endl;

	//calculate new grid info
	double h=mc._boun_box[1][0]-mc._boun_box[0][0];
	mc._boun_box[1][0]+=(new_size_x-mc._x_size)*h;
	mc._boun_box[1][1]+=(new_size_y-mc._y_size)*h;
	mc._boun_box[1][2]+=(new_size_z-mc._z_size)*h;
	mc._x_size=new_size_x;
	mc._y_size=new_size_y;
	mc._z_size=new_size_z;
	//cout<<mc._x_size<<" "<<mc._y_size<<" "<<mc._z_size<<endl;
	//cout<<"Here:::::"<<endl;
	//cout<<grid[102432]<<" "<<sign(grid[102432])<<endl;
	//cout<<grid[102433]<<" "<<sign(grid[102433])<<endl;
	//calculate intersection info 
	//x direction search
	int sss=0;
	for (int i = 0; i < mc._y_size; i++)
		for (int j = 0; j < mc._z_size; j++)
			for (int k = 0; k < mc._x_size - 1; k++)
			{
			  //if(i==68&&j==9&&k==24)
			  //{
			  //cout<<"asdasdadada"<<grid[i*mc._x_size+j*mc._x_size*mc._y_size+k]<<endl;
			  //}
			  if (sign(grid[i*mc._z_size + k*mc._z_size*mc._y_size + j]) != sign(grid[i*mc._z_size + (k+1)*mc._z_size*mc._y_size + j]))
				{
					std::vector<int> t_inte_dual;
					std::vector<double> t_inte_vert;
					t_inte_dual.push_back(k);
					t_inte_dual.push_back(i);
					t_inte_dual.push_back(j);
					t_inte_dual.push_back(k + 1);
					t_inte_dual.push_back(i);
					t_inte_dual.push_back(j);
					double t_x, t_y, t_z;
					t_x = (mc._boun_box[1][0] - mc._boun_box[0][0]) / (mc._x_size - 1)*(k + fabs(grid[i*mc._z_size + k*mc._z_size*mc._y_size + j]) / (fabs(grid[i*mc._z_size + k*mc._z_size*mc._y_size + j]) + fabs(grid[i*mc._z_size + (k+1)*mc._z_size*mc._y_size + j])));
					t_y = (mc._boun_box[1][1] - mc._boun_box[0][1]) / (mc._y_size - 1)*i;
					t_z = (mc._boun_box[1][2] - mc._boun_box[0][2]) / (mc._z_size - 1)*j;
					//cout<<t_x<<endl;
					t_inte_vert.push_back(t_x);
					t_inte_vert.push_back(t_y);
					t_inte_vert.push_back(t_z);
					mc._inte_dual.push_back(t_inte_dual);
					mc._inte_verts.push_back(t_inte_vert);
					//if(sss<1){
					  //cout<<grid[i*mc._x_size+j*mc._x_size*mc._y_size+k]<<endl;
					  //cout<<grid[i*mc._x_size+j*mc._x_size*mc._y_size+k+1]<<endl;
					//cout<<t_inte_dual[0]<<" "<<t_inte_dual[1]<<" "<<t_inte_dual[2]<<" "<<t_inte_dual[3]<<" "<<t_inte_dual[4]<<" "<<t_inte_dual[5]<<endl;
					//cout<<t_inte_vert[0]<<" "<<t_inte_vert[1]<<" "<<t_inte_vert[2]<<endl;sss++;}
				}
			}

	//y direction search
	for (int i = 0; i < mc._x_size; i++)
		for (int j = 0; j < mc._z_size; j++)
			for (int k = 0; k < mc._y_size - 1; k++)
			{
				if (sign(grid[j + i*mc._z_size*mc._y_size + k*mc._z_size]) != sign(grid[j + i*mc._z_size*mc._y_size + (k + 1)*mc._z_size]))
				{
					std::vector<int> t_inte_dual;
					std::vector<double> t_inte_vert;
					t_inte_dual.push_back(i);
					t_inte_dual.push_back(k);
					t_inte_dual.push_back(j);
					t_inte_dual.push_back(i);
					t_inte_dual.push_back(k + 1);
					t_inte_dual.push_back(j);
					double t_x, t_y, t_z;
					t_x = (mc._boun_box[1][0] - mc._boun_box[0][0]) / (mc._x_size - 1)*i;
					t_y = (mc._boun_box[1][1] - mc._boun_box[0][1]) / (mc._y_size - 1)*(k + fabs(grid[j + i*mc._z_size*mc._y_size + k*mc._z_size]) / (fabs(grid[j + i*mc._z_size*mc._y_size + k*mc._z_size]) + fabs(grid[j + i*mc._z_size*mc._y_size + (k + 1)*mc._z_size])));
					t_z = (mc._boun_box[1][2] - mc._boun_box[0][2]) / (mc._z_size - 1)*j;
					//cout<<t_x<<endl;
					t_inte_vert.push_back(t_x);
					t_inte_vert.push_back(t_y);
					t_inte_vert.push_back(t_z);
					mc._inte_dual.push_back(t_inte_dual);
					mc._inte_verts.push_back(t_inte_vert);
				}
			}

	//z direction search
	for (int i = 0; i < mc._x_size; i++)
		for (int j = 0; j < mc._y_size; j++)
			for (int k = 0; k < mc._z_size - 1; k++)
			{
				if (sign(grid[k + j*mc._z_size + i*mc._z_size*mc._y_size]) != sign(grid[k+1 + j*mc._z_size + i*mc._z_size*mc._y_size]))
				{
					std::vector<int> t_inte_dual;
					std::vector<double> t_inte_vert;
					t_inte_dual.push_back(i);
					t_inte_dual.push_back(j);
					t_inte_dual.push_back(k);
					t_inte_dual.push_back(i);
					t_inte_dual.push_back(j);
					t_inte_dual.push_back(k + 1);
					double t_x, t_y, t_z;
					t_x = (mc._boun_box[1][0] - mc._boun_box[0][0]) / (mc._x_size - 1)*i;
					t_y = (mc._boun_box[1][1] - mc._boun_box[0][1]) / (mc._y_size - 1)*j;
					t_z = (mc._boun_box[1][2] - mc._boun_box[0][2]) / (mc._z_size - 1)*(k + fabs(grid[k + j*mc._z_size + i*mc._z_size*mc._y_size]) / (fabs(grid[k + j*mc._z_size + i*mc._z_size*mc._y_size]) + fabs(grid[k+1 + j*mc._z_size + i*mc._z_size*mc._y_size])));
					//cout<<t_x<<endl;
					t_inte_vert.push_back(t_x);
					t_inte_vert.push_back(t_y);
					t_inte_vert.push_back(t_z);
					mc._inte_dual.push_back(t_inte_dual);
					mc._inte_verts.push_back(t_inte_vert);
				}
			}

	//grid info
	for (int i = 0; i < mc._z_size; i++)
		for (int j = 0; j < mc._y_size; j++)
			for (int k = 0; k < mc._x_size; k++)
			{
				std::vector<double> t_grid_info;
				t_grid_info.push_back(i);
				t_grid_info.push_back(j);
				t_grid_info.push_back(k);
				t_grid_info.push_back(grid[i + j*mc._z_size + k*mc._z_size*mc._y_size]);
				mc._grid_info.push_back(t_grid_info);
			}
	//initialize normal
	vector<double> init_normal;
	init_normal.push_back(0);
	init_normal.push_back(0);
	init_normal.push_back(0);
	for(int i=0; i<mc._inte_dual.size(); i++)
	  mc._inte_normal.push_back(init_normal);
	
	return 1;
}

int calculate_normal(MarchingCubes& mc)
{
  std::vector<int> divide;
  divide.assign(mc._vertices.size(),0);
  for(int i=0; i<mc._triangles.size(); i++)
    {
      //v1
      mc._inte_normal[mc._triangles[i].v1][0] += mc._triangle_normal[i][0];
      mc._inte_normal[mc._triangles[i].v1][1] += mc._triangle_normal[i][1]; 
      mc._inte_normal[mc._triangles[i].v1][2] += mc._triangle_normal[i][2];
      divide[mc._triangles[i].v1] += 1;
      
      //v2
      mc._inte_normal[mc._triangles[i].v2][0] += mc._triangle_normal[i][0]; 
      mc._inte_normal[mc._triangles[i].v2][1] += mc._triangle_normal[i][1]; 
      mc._inte_normal[mc._triangles[i].v2][2] += mc._triangle_normal[i][2];
      divide[mc._triangles[i].v2] += 1;

      //v3
      mc._inte_normal[mc._triangles[i].v3][0] += mc._triangle_normal[i][0];
      mc._inte_normal[mc._triangles[i].v3][1] += mc._triangle_normal[i][1];
      mc._inte_normal[mc._triangles[i].v3][2] += mc._triangle_normal[i][2];
      divide[mc._triangles[i].v3] += 1; 
    }

  for(int i=0; i<mc._vertices.size(); i++)
    {
      mc._vertices[i].nx = mc._inte_normal[i][0] / divide[i];
      mc._vertices[i].ny = mc._inte_normal[i][1] / divide[i];
      mc._vertices[i].nz = mc._inte_normal[i][2] / divide[i];
    }
  
  
  /*	for (int i = 0; i < mc._vertices.size(); i++)
	{
		std::vector<double> t_normal;
		t_normal.push_back(0);
		t_normal.push_back(0);
		t_normal.push_back(0);
		int count = 0;
		for (int j = 0; j < mc._triangles.size(); j++)
			if (i == mc._triangles[j].v1 || i == mc._triangles[j].v3 || i == mc._triangles[j].v3)
			{
				t_normal[0] += mc._triangle_normal[j][0];
				t_normal[1] += mc._triangle_normal[j][1];
				t_normal[2] += mc._triangle_normal[j][2];
				count++;
			}
		t_normal[0] /= count;
		t_normal[1] /= count;
		t_normal[2] /= count;
		mc._inte_normal.push_back(t_normal);
	}
  */
	return 1;
}
