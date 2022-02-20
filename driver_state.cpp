#include "driver_state.h"
#include <cstring>
using namespace std;

driver_state::driver_state()
{
}

driver_state::~driver_state()
{
    delete [] image_color;
    delete [] image_depth;
}

// This function should allocate and initialize the arrays that store color and
// depth.  This is not done during the constructor since the width and height
// are not known when this class is constructed.
void initialize_render(driver_state& state, int width, int height)
{
    state.image_width=width;
    state.image_height=height;
    //the image_color and image_depth are arrays, so fill out the arrays
    state.image_color= new pixel[width*height];
    state.image_depth= new float[width*height];

    for(int i = 0; i < (width*height); i++){

        //this initializes all pixels in image_color to black
        state.image_color[i] = make_pixel(0,0,0);
        //we can set image_depth to anything since we're ignoring it for now lab #5
        state.image_depth[i] = 100;
    }

    std::cout<<"TODO: allocate and initialize state.image_color and state.image_depth."<<std::endl;
}

// This function will be called to render the data that has been stored in this class.
// Valid values of type are:
//   render_type::triangle - Each group of three vertices corresponds to a triangle.
//   render_type::indexed -  Each group of three indices in index_data corresponds
//                           to a triangle.  These numbers are indices into vertex_data.
//   render_type::fan -      The vertices are to be interpreted as a triangle fan.
//   render_type::strip -    The vertices are to be interpreted as a triangle strip.
void render(driver_state& state, render_type type)
{
    std::cout<<"TODO: implement rendering."<<std::endl;
    
    switch(type){
        case render_type::triangle:{
            
            //first step: allocal an array of data_geometry objects, one for each vertex
            data_geometry objs[3];//there is 3 vertices per triangle
            data_vertex vertex[3]; // 3 vertices per triangle, so 1 data_vertex per vertice
            for (int i =0 ; i<3 ; ++i){ //we're going to go through each vertex in a loop and allocate the data member float* data
                //for each vertex i, we're going to fill in for each of the j floats in each vertex
                //j is incremented by 3*state.floats_per_vertex because we are going by vertex
                for (int j = 0 ; j< state.num_vertices*state.floats_per_vertex; j += 3 *state.floats_per_vertex){
                    vertex[i].data= &state.vertex_data[(state.floats_per_vertex * i) + j];
                }
                //we need to fill in the contents of the objects, using vertex shader
                state.vertex_shader(vertex[i], objs[i], state.uniform_data);
                //we now have on data_gepmetry object for each vertex fully filled
            }
            //call rasterize_triangle to rasterize the triangle
            rasterize_triangle(state, objs[0], objs[1], objs[2]); 
        }
        break;
        //leave the remaining cases empty because we don't need them for now
        case render_type::indexed:{}
        break;
        case render_type::strip: {}
        break;
        case render_type::fan: {}
        break;
    }
}


// This function clips a triangle (defined by the three vertices in the "in" array).
// It will be called recursively, once for each clipping face (face=0, 1, ..., 5) to
// clip against each of the clipping faces in turn.  When face=6, cli triangle should
// simply pass the call on to rasterize_triangle.
void clip_triangle(driver_state& state, const data_geometry& v0,
    const data_geometry& v1, const data_geometry& v2,int face)
{
    if(face==6)
    {
        rasterize_triangle(state, v0, v1, v2);
        return;
    }
    std::cout<<"TODO: implement clipping. (The current code passes the triangle through without clipping them.)"<<std::endl;
    clip_triangle(state, v0, v1, v2,face+1);
}

// Rasterize the triangle defined by the three vertices in the "in" array.  This
// function is responsible for rasterization, interpolation of data to
// fragments, calling the fragment shader, and z-buffering.
void rasterize_triangle(driver_state& state, const data_geometry& v0,
    const data_geometry& v1, const data_geometry& v2)
{
    std::cout<<"TODO: implement rasterization"<<std::endl;

    //first we fill out the coordinates for each vertex
    float A_x = ((state.image_width/2) * v0.gl_Position[0]);
    float A_y = ((state.image_height/2) * v0.gl_Position[1]);
    float B_x = ((state.image_width/2) * v1.gl_Position[0]);
    float B_y = ((state.image_height/2) * v1.gl_Position[1]);
    float C_x = ((state.image_width/2) * v2.gl_Position[0]);
    float C_y = ((state.image_height/2) * v2.gl_Position[1]);
    //these are the coordinates for point P where i is x and j is y
    float  P_i = (state.image_width/2) - 0.5;
    float  P_j = (state.image_height/2) - 0.5;
    //These are the coordinates that relate the vertices to the point P coordinates
    //we need to transfer the data_geometry positions from NDC's to pixel coordinates
    float Ax_i = A_x + P_i;
    float Ay_j = A_y + P_j;
    float Bx_i = B_x + P_i;
    float By_j = B_y + P_j;
    float Cx_i = C_x + P_i;
    float Cy_j = C_y + P_j;
    //this was a given formula from the lab sheet
    float AREA_ABC = 0.5 * (((Bx_i * Cy_j - (Cx_i * By_j)) - ((Ax_i* Cy_j) - (Cx_i * Ay_j)) + ((Ax_i* By_j) - (Bx_i * Ay_j))));

    //we have a forloop for triangle rasterization for all x from xmin to xmax, with a forloop for all y from ymin to ymax
    //i put using namespace std at the top so the max and min functions would work
    for (int i = min(min(Ax_i,Bx_i),Cx_i); i < max(max(Ax_i,Bx_i),Cx_i); i++){
        for (int j=min(min(Ay_j,By_j),Cy_j); j<max(max(Ay_j,By_j),Cy_j) ; j++){
            //for alpha beta and gamma you just use the same equation for ABC but swap out the letter youre solving for with P
            //or in this case, the i and j's which are going through each min to max values
            float alpha = 0.5 * (((Bx_i * Cy_j)-(Cx_i * By_j))-(( i * Cy_j) - (Cx_i * j))+(( i * By_j)-(Bx_i * j))) /AREA_ABC;
            float beta = 0.5 * ((( i * Cy_j)-(Cx_i * j))-((Ax_i * Cy_j) - (Cx_i * Ay_j))+((Ax_i * j)-( i * Ay_j))) /AREA_ABC;
            float gamma = 0.5 * (((Bx_i * j)-( i * By_j))-((Ax_i * j)-( i * Ay_j))+((Ax_i * By_j) - (Bx_i * Ay_j))) /AREA_ABC;
        //this implies 
        if (alpha >=0 && beta >=0 && gamma <=1){
            int index = (state.image_width  * j) + i;
            state.image_color[index] = make_pixel(255,255,255);
        }
        }
   }

}

