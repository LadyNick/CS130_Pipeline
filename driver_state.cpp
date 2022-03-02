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
        state.image_depth[i] = 100.0;
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
            data_geometry objs[3];//there is 3 vertices per triangle, but we'll just use vertices per tria
            data_vertex vertex[3]; // 3 vertices per triangle, so 1 data_vertex per vertice
            for (int i = 0 ; i<state.num_vertices*state.floats_per_vertex; i += 3*state.floats_per_vertex){ 
                for (int j = 0 ; j< 3; j++){
                    //fill in the data for each vertices
                    vertex[j].data= &state.vertex_data[(state.floats_per_vertex * j) + i];
                    objs[j].data=vertex[j].data;
                    //we need to fill in the contents of the objects, using vertex shader
                    state.vertex_shader(vertex[j], objs[j], state.uniform_data);
                }
                //we now have on data_gepmetry object for each vertex fully filled
            //call rasterize_triangle to rasterize the triangle
            rasterize_triangle(state, objs[0], objs[1], objs[2]); 
            }
        }
        break;
        //leave the remaining cases empty because we don't need them for now
        case render_type::indexed:{}
        break;
        case render_type::strip: {}
        break;
        case render_type::fan: {}
        break;
        case render_type::invalid:{}
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
    //these are the coordinates for pixel P where i is x and j is ... because its in the center of every pixel u could say
    //so its like if we had a grid, normally 0,0 would be 0,0 the origin but we want 0,0 to become 0.5 0.5
    float  width_alter = (state.image_width/2) - 0.5;
    float  height_alter = (state.image_height/2) - 0.5;
    //These are the coordinates that relate the vertices to the pixel coordinates
    //we need to transfer the data_geometry positions from NDC's to pixel coordinates
    //so theyre in the center of every pixel
    float Ax = A_x + width_alter;
    float Ay = A_y + height_alter;
    float Bx = B_x + width_alter;
    float By = B_y + height_alter;
    float Cx = C_x + width_alter;
    float Cy = C_y + height_alter;
    //this sets up the bounding box
    int minx = min(min(Ax,Bx),Cx);
    int maxx = max(max(Ax,Bx),Cx);
    int miny = min(min(Ay,By),Cy);
    int maxy = max(max(Ay,By),Cy);
    //this was a given formula from the lab sheet
    float AREA_ABC = 0.5 * ((Bx * Cy - Cx * By) + (Cx * Ay - Ax * Cy) + (Ax * By - Bx * Ay));

    //we have a forloop for triangle rasterization for all x from xmin to xmax, with a forloop for all y from ymin to ymax
    //i put using namespace std at the top so the max and min functions would work
    for (int i = minx; i < maxx; i++){ //track of current x
        for (int j=miny; j < maxy ; j++){ //track of current y
            //for alpha beta and gamma you just use the same equation for ABC but swap out the letter youre solving for with i or p
            //or in this case, the i and j's which are going through the bounding box
            float alpha = 0.5 * ((Bx * Cy - Cx * By) + (Cx * j - i * Cy) + (i * By - Bx * j)) /AREA_ABC;
            float beta = 0.5 * (( i * Cy - Cx * j) + (Cx * Ay - Ax * Cy) + (Ax * j - i * Ay)) /AREA_ABC;
            float gamma = 1.0 - alpha - beta;
            //if alpha and beta are positive then gamma must be positive 
        if(alpha >= 0 && beta >= 0 && gamma >= 0){


            //if the pixel is inside the triangle, then do interpolation
            data_fragment fragment;
            float frag_data[state.floats_per_vertex];
            fragment.data = frag_data;

            for(int k = 0; k < state.floats_per_vertex; k++){

                switch(state.interp_rules[k]){
                     case interp_type::invalid:{
                        //invalid interpolation
                    }
                    break;
                    case interp_type::flat:{
                        fragment.data[k]= v0.data[k];
                    }
                    break;
                    case interp_type::smooth:{
                       //smooth interpolation    
                    }
                    break;
                    case interp_type::noperspective:{   
                        fragment.data[k] = alpha * v0.data[k] + beta * v1.data[k]+gamma * v2.data[k];
                    }
                    break;
                }
            }
        //call fragment shader for the fragment we just interpolated
        data_output out;
        state.fragment_shader(fragment, out, state.uniform_data);

        //z buffer, the ones that are closer to the front are the ones we want to see on the image
        float z_buff = (alpha * v0.gl_Position[2]/v0.gl_Position[3])+(beta*v1.gl_Position[2]/v1.gl_Position[3])+(gamma*v2.gl_Position[2]/v2.gl_Position[3]);
        int index = (state.image_width * j) + i;

        //we will not color the indexed pixel if it doesn't meet the requirement
        if(state.image_depth[index] > z_buff){
          state.image_depth[index] = z_buff;
          //establish color for the pixel that is the frontmost 
          state.image_color[index] = make_pixel(255 * out.output_color[0], 255* out.output_color[1], 255* out.output_color[2]);
        }
        }
        }
   }

}

