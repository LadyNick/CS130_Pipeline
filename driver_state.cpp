#include "driver_state.h"
#include <cstring>
using namespace std;

void makevertex(data_geometry& vertexmake, const data_geometry& verta, const data_geometry& vertb, int face, driver_state& state, float* new_data);

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
                    objs[j].data = vertex[j].data= &state.vertex_data[(state.floats_per_vertex * j) + i];
                    //we need to fill in the contents of the objects, using vertex shader
                    state.vertex_shader(vertex[j], objs[j], state.uniform_data);
                }
                //we now have on data_gepmetry object for each vertex fully filled
            //call rasterize_triangle to rasterize the triangle
            //rasterize_triangle(state, objs[0], objs[1], objs[2]); 
            clip_triangle(state, objs[0], objs[1], objs[2], 0);
            }
        }
        break;
        //leave the remaining cases empty because we don't need them for now
        case render_type::indexed:{
            data_geometry objs[3];
            data_vertex vertex[3];
            for(int i = 0; i<3*state.num_triangles; i+=3){
                for(int j=0; j<3; j++){
                    objs[j].data = vertex[j].data=&state.vertex_data[state.index_data[i+j]*state.floats_per_vertex];
                    state.vertex_shader(vertex[j], objs[j], state.uniform_data);
                }
                //rasterize_triangle(state, objs[0], objs[1], objs[2]);
                //clipping triangle goes here
                clip_triangle(state, objs[0], objs[1], objs[2], 0);
            }    
        }
        break;
        case render_type::strip: {
            data_geometry objs[3];
            data_vertex vertex[3];
            for(int i = 0; i < (state.num_vertices-2); i++){
                for(int j =0; j<3; j++){
                    objs[j].data = vertex[j].data = &state.vertex_data[state.floats_per_vertex * (i+j)];
                    state.vertex_shader(vertex[j],objs[j],state.uniform_data);
                }
                //rasterize_triangle(state, objs[0], objs[1], objs[2]);
                //clipping goes here
                clip_triangle(state, objs[0], objs[1], objs[2], 0);
            }
        }
        break;
        case render_type::fan: {
            data_geometry objs[3];
            data_vertex vertex[3];
            for(int i = 0; i<state.num_vertices; i++){
                for(int j = 0; j<3; j++){
                    int k = ((i+j)*state.floats_per_vertex );
                    if(j == 0){
                        k=0;
                    }
                    objs[j].data = vertex[j].data = state.vertex_data + k;
                    state.vertex_shader(vertex[j], objs[j], state.uniform_data);
                    }
                }
                //rasterize_triangle(state, objs[0], objs[1], objs[2]);
                //clipping
                clip_triangle(state, objs[0], objs[1], objs[2], 0);
            }
        break;
        case render_type::invalid: {
            cout << "Invalid" << endl; //i think we had an invalid case but i dunno if 
            //we were supposed to do anything specific for it, i guess for testing perhaps
        }
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
    data_geometry arr[3];
  data_geometry o[3];
  data_vertex ver[3];
    data_geometry g[3];
 
        arr[0].gl_Position= v0->gl_Position;
        arr[1].gl_Position = v1->gl_Position;
        arr[2].gl_Position = v2->gl_Position;
     
  int  position =0, value =0;

    if(face==6)
    {

        rasterize_triangle(state, v0, v1, v2);
        return;
    }
      else if (face==0 )
    {
      position = 0;
      value =1;
    }
      else if (face ==1)// set the position values
    {
        position =0;
        value =-1;
    }
    else if (face == 2)
    {
      position =1;
      value =1;
    }
    else if (face == 3 )
    {
      position =1;
      value =-1;
    }
    else if (face == 4 )
    {
      position =2;
      value =1;
    }
    else if (face == 5 )
    {
      position =2;
       value =-1;
    }
   
vec3 tri = {0,0,0};
int count =0;

  if (value == -1)
    {
      if(v0->gl_Position[position] >= -v0->gl_Position[3])
        {
        count++;
        tri[0] = 1;
        }
     }
    else if(v0->gl_Position[position] <=  v0->gl_Position[3])
     {
      count++;
      tri[0] = 1;
      }
    if (value == -1)
    {
      if(v1->gl_Position[position] >= -v1->gl_Position[3])
        {
        count++;
        tri[1] = 1;
        }
     }
    else if(v1->gl_Position[position] <=  v1->gl_Position[3])
     {
      count++;
      tri[1] = 1;
      }
    if (value == -1)
    {
      if(v2->gl_Position[position] >= -v2->gl_Position[3])
        {
        count++;
        tri[2] = 1;
        }
     }
    else if(v2->gl_Position[position] <=  v2->gl_Position[3])
     {
      count++;
      tri[2] = 1;
      }
    
    
if (count ==3 )
{
//everything inside
  clip_triangle(state,v0,v1,v2 ,face+1);
  return;
}
if (count ==0 )
{
  //everything is outside, ..... so no cares
  return;
}

data_geometry* in[3];
in[0] = v0;
in[1] = v1;
in[2] = v2;
    
const data_geometry a =0;
const data_geometry b =0;
const data_geometry c = 0;
for (int i =0; i<3; i++)
{
  if ((count ==1 && tri[i]==1 ) || (count ==2 && tri[i]==0))
   {
//setting the triangles vertices;
     a = in[i];
     b = in[(i+1)%3];
     c = in[(i+2)%3];
   }
}
float position_A = a->gl_Position[3];
float position_B = b->gl_Position[3];
float position_C = c->gl_Position[3];
float beta = ((value * position_B) - b->gl_Position[position]) / ((a->gl_Position[position] - (value * position_A)) + ((value * position_B) - b->gl_Position[position]));
float gamma = ((value * position_C) - c->gl_Position[position]) / ((a->gl_Position[position] - (value * position_A)) + ((value * position_C) - c->gl_Position[position]));

data_geometry* ba = new data_geometry();
data_geometry* ca = new data_geometry();
 float BA_data[state.floats_per_vertex];
ba->data= BA_data;
float CA_data[state.floats_per_vertex];
ca->data= CA_data;
for (int j =0; j<state.floats_per_vertex;  j++)
{
  switch(state.interp_rules[j])
  {
    case interp_type::flat:
    {
      ba->data[j]=a->data[j];
      ca->data[j]=a->data[j];
    }break;
    case interp_type::smooth:
    {
      ba->data[j]= beta* a->data[j]+ ((1-beta)*b->data[j]);
      ca->data[j]= gamma* a->data[j]+((1-gamma)*c->data[j]);

    }break;
    case interp_type::noperspective:
    {
      float  k = (beta*position_A)+((1-beta)*position_B);
      float alpha_prime = (beta * position_A)/k;
      ba->data[j]= (alpha_prime * a->data[j])+((1-alpha_prime )*b->data[j]);
      float k_gamma = (gamma *position_A)+((1-gamma)*position_C);
      float gamma_prime = (gamma*position_A)/k_gamma;
      ca->data[j]= (gamma_prime*a->data[j])+((1-gamma_prime)*c->data[j]);
    }break;
    default : std::cout<<"GOD!HELP ME!!!!!!"<<std::endl;

  }

}
for (int i = 0;i<4; i++)
{
  ba->gl_Position[i]= (beta * a->gl_Position[i]) + ((1-beta)*b->gl_Position[i]);
  ca->gl_Position[i]= (gamma * a->gl_Position[i]) + ((1-gamma)*c->gl_Position[i]);

}


if (count ==1)
{
  g[0]=a;
  g[1]=ba;
  g[2]=ca;
  clip_triangle(state,g[0], g[1], g[2],face+1);

  return;
}
if (count ==2 )
{
  g[0]=ba;
  g[1]=b;
  g[2]=c;
  clip_triangle(state,g[0], g[1], g[2],face+1);
  g[0]=c;
  g[1]=ca;
  g[2]=ba;
  clip_triangle(state,g[0], g[1], g[2],face+1);
  return;
    
    

   
   /* if(face==6)
    {//this means that this is the default case
        rasterize_triangle(state, v0, v1, v2);
        //i might try removing this part since technically there is 5 faces, 
        //not 6 and we're starting from 0
        return;
    }
    //std::cout<<"TODO: implement clipping. (The current code passes the triangle through without clipping them.)"<<std::endl;
    //clip_triangle(state, v0, v1, v2,face+1);    this will come later

    //first off there are 6 faces, 2 for +-x, 2 for +-y, 2 for +-x, so 01, 23, 45
    //make an if for each case to handle the plane
    
    //then for each plane we have to determine what vertices are in or out
    //then depending on which ones are in or out we have 4 cases

    //all inside, pass triangle to next stage
    //all outside, do nothing, discard the triangle (which i assume means return)
    //1 in 2 out --> find intersection points, then do a new triangle 
            //with the one that is in, and the 2 intersecting with the plane as the edge
            //then pass that triangle to the next clipping stage
    //2 in 1 out --> find the intersection points, then  make 2 triangles, each differing by 1
            //a ab bc      a bc c
    

    bool avert = 0;
    bool bvert = 0;
    bool cvert = 0;

    //this will be used to determine which vertices are inside or outside

    //negative x plane
    if(face == 0){ 
        //inside if -w <= x
        //the w value in in gl_position 3
        if(v0.gl_Position[0] >= -1 * v0.gl_Position[3]){
            avert = true;
        }
        else{ avert = false; }
        if(v1.gl_Position[0] >= -1 * v1.gl_Position[3]){
            bvert = true;
        } 
        else{ bvert = false; }
        if(v2.gl_Position[0] >= -1 * v2.gl_Position[3]){
            cvert = true;
        }
        else{ cvert = false; }
    }
    //positive x plane
    if(face == 1){
        //inside if x <= w
        //the w value in in gl_position 3
        if(v0.gl_Position[0] <= v0.gl_Position[3]){
            avert = true;
        }
        else{ avert = false; }
        if(v1.gl_Position[0] <= v1.gl_Position[3]){
            bvert = true;
        } 
        else{ bvert = false; }
        if(v2.gl_Position[0] <= v2.gl_Position[3]){
            cvert = true;
        }
        else{ cvert = false; }
    }
    //negative y plane
    if(face == 2){
        //inside if -w <= y
        //the w value in in gl_position 3
        if(v0.gl_Position[1] >= -1 * v0.gl_Position[3]){
            avert = true;
        }
        else{ avert = false; }
        if(v1.gl_Position[1] >= -1 * v1.gl_Position[3]){
            bvert = true;
        } 
        else{ bvert = false; }
        if(v2.gl_Position[1] >= -1 * v2.gl_Position[3]){
            cvert = true;
        }
        else{ cvert = false; }
    }
    //positive y plane
    if(face == 3){
        //inside if y <= w
        //the w value in in gl_position 3
        if(v0.gl_Position[1] <= v0.gl_Position[3]){
            avert = true;
        }
        else{ avert = false; }
        if(v1.gl_Position[1] <= v1.gl_Position[3]){
            bvert = true;
        } 
        else{ bvert = false; }
        if(v2.gl_Position[1] <= v2.gl_Position[3]){
            cvert = true;
        }
        else{ cvert = false; }
    }
    //negative z plane
    if(face == 4){
        cout << "face4" << endl;
        cout << v0.gl_Position[0] << v0.gl_Position[1] << v0.gl_Position[2] << v0.gl_Position[3] << endl;
        cout << v1.gl_Position[0] << " " << v1.gl_Position[1] <<" "<< v1.gl_Position[2] << endl;
        cout << v2.gl_Position[0] <<" " << v2.gl_Position[1] << " " << v2.gl_Position[2] << endl;

        //inside if -w <= z
        //the w value in in gl_position 3
        if(v0.gl_Position[2] >= -1 * v0.gl_Position[3]){
            avert = true;
        }
        else{ avert = false; }
        if(v1.gl_Position[2] >= -1 * v1.gl_Position[3]){
            bvert = true;
        } 
        else{ bvert = false; }
        if(v2.gl_Position[2] >= -1 * v2.gl_Position[3]){
            cvert = true;
        }
        else{ cvert = false; }
        cout << avert << bvert << cvert << endl;
    }
   
    //positive z plane
    if(face == 5){
      cout << "got to face 5" << endl;
        //inside if z <= w
        //the w value in in gl_position 3
        if(v0.gl_Position[2] <= v0.gl_Position[3]){
            avert = true;
        }
        else{ avert = false; }
        if(v1.gl_Position[2] <= v1.gl_Position[3]){
            bvert = true;
        } 
        else{ bvert = false; }
        if(v2.gl_Position[2] <= v2.gl_Position[3]){
            cvert = true;
        }
        else{ cvert = false; }
        cout << avert << bvert << cvert << endl;
    }

    //at this point we know what vertices are inside and what vertices are outside
    //meaning now we should make triangles if there are any triangles to be made
    //HOWEVER, we don't know how many triangles will be made and we don't know 
    //what vertices will make the edges if there are any triangles, so we have 
    //to try every options
    data_geometry tri1[3];
    data_geometry tri2[3];
    float data1[MAX_FLOATS_PER_VERTEX];
    float data2[MAX_FLOATS_PER_VERTEX];
    //if all are inside, pass the triangle to the next stage
    if(avert && bvert && cvert){
        clip_triangle(state, v0, v1, v2, face+1);
        
    }
    //if all are outside, discard the triangle
    else if(!avert && !bvert && !cvert){
        //do nothing? return maybe
        return;
    }
    //now we need to handle the remaning cases that result in triangles
    //2 out 1 in, 1 triangle
    else if(!avert && !bvert && cvert){
        //lets make a triangle
        //c is base because c is only one in
        //data_geometry tri1[3];
        tri1[0].gl_Position = v2.gl_Position;
        tri1[0].data = v2.data;

        //and now the annoying part where we have to make the other vertices
        //ca and cb, which is gonna be a while so I'll make a helper so I don't
        //have to rewrite it, since theres more interpolation before clipping again
        makevertex(tri1[1], v2, v0, face, state, data1);
        makevertex(tri1[2], v2, v1, face, state, data2);
        //since at this moment the triangle is changeable, in order to
        //pass it through clip again we need to mkae it a constant
        //const data_geometry trigl1 = const_cast<data_geometry>(tri1);
        clip_triangle(state, tri1[0], tri1[1], tri1[2], face+1);

    }
    else if(!avert && bvert && !cvert){
        //lets make a triangle
        //use b as the base since b is in
        //data_geometry tri1[3];
        tri1[0].gl_Position = v1.gl_Position;
        tri1[0].data = v1.data;
        //make vertices ba and bc
        makevertex(tri1[1], v1, v2, face, state, data1);
        makevertex(tri1[2], v1, v0, face, state, data2);
        //since at this moment the triangle is changeable, in order to
        //pass it through clip again we need to mkae it a constant
        //const data_geometry trigl1 = const_cast<data_geometry>(tri1);
        clip_triangle(state, tri1[0], tri1[1], tri1[2], face+1);

    }
    else if(avert && !bvert && !cvert){
        //lets make a triangle
        //since a is the one in, we're working with a
       // data_geometry tri1[3];
        tri1[0].gl_Position = v0.gl_Position;
        tri1[0].data = v0.data;
        //make vertices ab and ac
        makevertex(tri1[1], v0, v1, face, state, data1);
        makevertex(tri1[2], v0, v2, face, state, data2);
        //since at this moment the triangle is changeable, in order to
        //pass it through clip again we need to mkae it a constant
        //const data_geometry trigl1 = const_cast<data_geometry>(tri1);
        clip_triangle(state, tri1[0], tri1[1], tri1[2], face+1);
    }
    //now we can handle the 1out 2 in cases
    else if(avert && bvert && !cvert){
        //data_geometry tri1[3];
        //data_geometry tri2[3];
        //this time because we have 2 in, we can set two of the vertices to the ones that are in
        tri1[0].gl_Position=v0.gl_Position; //v0 is for a
        tri1[0].data = v0.data;
        tri1[1].gl_Position=v1.gl_Position;
        tri1[1].data=v1.data;  //v1 is for b
        makevertex(tri1[2],v0,v2,face,state, data1);
        tri2[0].gl_Position=v1.gl_Position; //v1 is for b
        tri2[0].data=v1.data;
        makevertex(tri2[1], v1, v2, face, state, data2);
        tri2[2].gl_Position = tri1[2].gl_Position;
        tri2[2].data = tri1[2].data;
        clip_triangle(state, tri1[0], tri1[1], tri1[2], face+1);
        clip_triangle(state, tri2[0], tri2[1], tri2[2], face+1);
    }
    else if(avert && !bvert && cvert){
        //data_geometry tri1[3];
        //data_geometry tri2[3];
        //this time because we have 2 in, we can set two of the vertices to the ones that are in
        tri1[0].gl_Position=v2.gl_Position; //v2 is for c
        tri1[0].data = v2.data;
        tri1[1].gl_Position=v0.gl_Position;
        tri1[1].data=v0.data;  //v0 is for a
        makevertex(tri1[2],v2,v1,face,state, data1);
        tri2[0].gl_Position=v0.gl_Position; //v0 is for a
        tri2[0].data=v0.data;
        makevertex(tri2[1], v0, v1, face, state, data2);
        tri2[2].gl_Position = tri1[2].gl_Position;
        tri2[2].data = tri1[2].data;
        clip_triangle(state, tri1[0], tri1[1], tri1[2], face+1);
        clip_triangle(state, tri2[0], tri2[1], tri2[2], face+1);
    }
    else if(!avert && bvert && cvert){
        //data_geometry tri1[3];
        //data_geometry tri2[3];
        //this time because we have 2 in, we can set two of the vertices to the ones that are in
        tri1[0].gl_Position=v1.gl_Position; //v1 is for b
        tri1[0].data = v1.data;
        tri1[1].gl_Position=v2.gl_Position;
        tri1[1].data=v2.data;  //v2 is for c
        makevertex(tri1[2],v2,v1,face,state, data1);
        tri2[0].gl_Position=v2.gl_Position; //v2 is for c
        tri2[0].data=v2.data;
        makevertex(tri2[1], v2, v0, face, state,data2);
        tri2[2].gl_Position = tri1[2].gl_Position;
        tri2[2].data = tri1[2].data;
        cout << tri1[2].gl_Position[0] << " " << tri1[2].gl_Position[1] <<" "<< tri1[2].gl_Position[2] << endl;
        cout << tri2[1].gl_Position[0] << " " << tri2[1].gl_Position[1] <<" "<< tri2[1].gl_Position[2] << endl;
        //cout << v1.gl_Position[0] << " " << v1.gl_Position[1] <<" "<< v1.gl_Position[2] << endl;
        clip_triangle(state, tri1[0], tri1[1], tri1[2], face+1);
        clip_triangle(state, tri2[0], tri2[1], tri2[2], face+1);
    }
    
}

void makevertex(data_geometry& vertexmake, const data_geometry& verta, const data_geometry& vertb, int face, driver_state& state, float* new_data){
    
    //from lecture, we got an equation for positive alpha and negative alpha
    //alpha left = -(xb + wb)/ ( (xa + wa) - (xb + wb) )
    //alpha right = (wb - xb) / ( (xa - wa) - (xb - wb) )
    //we'll have a vertex a and a vertex b just to go with the equation
    //so first we calculate alpha
    float alphapers;
    

    if(face == 1 || face == 3 || face == 5){//its from the right, so alpha right
        alphapers = (vertb.gl_Position[3] - vertb.gl_Position[face/2]) / ( (verta.gl_Position[face/2] - verta.gl_Position[3]) - (vertb.gl_Position[face/2] - vertb.gl_Position[3]) );
    }
    else{//face is soemthing else, so its negative
        alphapers = -(vertb.gl_Position[face/2] + vertb.gl_Position[3]) / (  (verta.gl_Position[face/2] + verta.gl_Position[3]) - (vertb.gl_Position[face/2] + vertb.gl_Position[3]) );
    }
    
    //from lecture, p = alpha*a + (1-alpha)*b
    //now we have the correct position in the vertex were making
    vertexmake.gl_Position = alphapers * verta.gl_Position + (1 - alphapers) * vertb.gl_Position;
    

    //mnow were gonna do a repeat of what we did in rasterize triangle where we did interpolation
    for(int i =0; i<state.floats_per_vertex; i++){
      
        switch(state.interp_rules[i]){
            case interp_type::flat:{
               
                new_data[i] = verta.data[i];
           

            }
            break;
            case interp_type::smooth:{
             
                new_data[i] = alphapers * verta.data[i] + (1-alphapers)*vertb.data[i];
                
            }
            break;
            case interp_type::noperspective:{
                //different alpha 
             
                float nopersa = alphapers * verta.gl_Position[3]/(alphapers * verta.gl_Position[3] + (1-alphapers)*vertb.gl_Position[3]);
                new_data[i] = nopersa*verta.data[i] + (1-nopersa) *vertb.data[i];

            }
            break;
            case interp_type::invalid:{
                //dont do anything?
            }
            break;
        }

    }
    vertexmake.data = new_data;*/
}

// Rasterize the triangle defined by the three vertices in the "in" array.  This
// function is responsible for rasterization, interpolation of data to
// fragments, calling the fragment shader, and z-buffering.
void rasterize_triangle(driver_state& state, const data_geometry& v0,
    const data_geometry& v1, const data_geometry& v2)
{
    //std::cout<<"TODO: implement rasterization"<<std::endl;

    //first we fill out the coordinates for each vertex
    float A_x = ((state.image_width/2) * v0.gl_Position[0])/v0.gl_Position[3];
    float A_y = ((state.image_height/2) * v0.gl_Position[1])/v0.gl_Position[3];
    float B_x = ((state.image_width/2) * v1.gl_Position[0])/v1.gl_Position[3];
    float B_y = ((state.image_height/2) * v1.gl_Position[1])/v1.gl_Position[3];
    float C_x = ((state.image_width/2) * v2.gl_Position[0])/v2.gl_Position[3];
    float C_y = ((state.image_height/2) * v2.gl_Position[1])/v2.gl_Position[3];
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
                        cout << "Invalid" << endl; //not too sure what to do here but I don't think it's too important
                    }
                    break;
                    case interp_type::flat:{
                        fragment.data[k]= v0.data[k];
                    }
                    break;
                    case interp_type::smooth:{
                       //smooth interpolation, has to be perspective correct
                       //we convert screen space barycentric coordiantes to world space coordinates
                       float m = alpha/v0.gl_Position[3] + beta/v1.gl_Position[3] + gamma/v2.gl_Position[3];
                       float smoothalpha = alpha/(v0.gl_Position[3] * m);
                       float smoothbeta = beta/(v1.gl_Position[3] * m);
                       float smoothgamma = gamma/(v2.gl_Position[3] * m);
                       //now interpolate
                       fragment.data[k] = smoothalpha * v0.data[k] + smoothbeta * v1.data[k] + smoothgamma * v2.data[k];  
                    }
                    break;
                    case interp_type::noperspective:{   
                        fragment.data[k] = alpha * v0.data[k] + beta * v1.data[k] + gamma * v2.data[k];
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

