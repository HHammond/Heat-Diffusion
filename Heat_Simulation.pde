/*//////////////////////////////////////////////////////////////////////////////

    Copyright (C) 2011  Henry Hammond
    email: HenryHHammond92@gmail.com
           HenryHammond@cmail.carleton.ca

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    For a copy of the GNU Lesser General Public License, see
    <http://www.gnu.org/licenses/>.
    
    /////////////////////////////////////////////////////////////////////////////
    
    This software was written for the 2013 MCM Problem A to simulate heat 
    distribution within differently shaped baking trays. This software assumes 
    the trays to be made from an alluminum alloy (as most kitchenware is) and
    makes the assumption that the tray is put into a uniform convection oven 
    filled with regular air. 

//////////////////////////////////////////////////////////////////////////////*/


import java.awt.Polygon;
import java.util.ArrayList;

int simWidth = 1280/4;
int yOffset = 20;

//meshing data
int w=simWidth;
int l=simWidth;
int h=8;

//meshing scale factor
int delx = 1;
int dely = 1;
int delz = 1;

//time sensitive data
float timestep = .1;
float time = 0;
//material data

//aluminum

//density
float row = 270.0; // kg/m^3

//specific heat capacity
float C = 0.902*1000; // J/kgK

//thermal conductivity
float k = 205; // W/mK

//heat transfer coefficient
float heatTransfer = 0.1;



//Oven data
float airTemp = 175;
float surfaceTemp = 20;

float maxTime = 30;

//graphing data
float maxTemp = airTemp+10;
float minTemp = 0;

//Other globals
boolean paused = false;

//simulations
Simulation sim1;
Simulation sim2;
Simulation sim3;
Simulation sim4;


PGraphics heatGrad;

void setup() {
  
  // Initialize window
  size(4*simWidth, 2*simWidth+yOffset);

  // set high framerate, this doesn't help performance much
  frameRate(120);
  
  // initialize number of sides per shape in simulation
  int n1 = 3;
  int n2 = 4;
  int n3 = 6;
  int n4 = 20;

  //create standard ratio to compare
  float w = 150;
  
  //set areas of each polygon
  float ratio12 = sqrt( (n2/n1)*(tan(PI/n1)/tan(PI/n2)) );
  float ratio13 = sqrt( (n3/n1)*(tan(PI/n1)/tan(PI/n3)) );
  float ratio14 = sqrt( (n4/n1)*(tan(PI/n1)/tan(PI/n4)) );
    
  // initialize simulations
  sim1 = new Simulation(n1,simWidth*0,yOffset,w);
  sim2 = new Simulation(n2,simWidth*1,yOffset,w/ratio12);
  sim3 = new Simulation(n3,simWidth*2,yOffset,w/ratio13);
  sim4 = new Simulation(n4,simWidth*3,yOffset,w/ratio14);

  /*
  //Cusom polygon if desired
  Polygon k = new Polygon();
  int rw = 70;
  k.addPoint(100,110);
  k.addPoint(100+rw,110);
  k.addPoint(100+rw,110+2*rw);
  k.addPoint(100,110+2*rw);
  
  sim1 = new Simulation(k,simWidth*0,yOffset);
  */

  
  //create graphics for heat gradient scale
  heatGrad = createGraphics(255,10);
  heatGrad.beginDraw();  
  heatGrad.noStroke();
  for(int i=0;i<256;i++){
    heatGrad.fill(i);
    heatGrad.rect(i,0,1,10);
  }
  heatGrad.endDraw();
  
}


void draw() {
    
  if(!paused){
    time+=timestep;
    background(200);
    
    //update simulations
    sim1.draw();
    sim2.draw();
    sim3.draw();
    sim4.draw();
    
    //update and draw graphcs
    int graphX = 0;
    int graphY = simWidth+yOffset;
    
    color c1 = color(255,0,0);    //  line 1 - inner temperature
    color c2 = color(0,0,255);    //  line 2 - outer temperature
    color c3 = color(255,127,0);  //  line 3 - temperature difference
    
    //divider between simulation screens
    color divider = color(150);
    
    //draw grpahs for each simulation
    //sim 1
    Graph g1 = new Graph(graphX,graphY,simWidth,simWidth,sim1.inner.values,c1,0,maxTemp);
    Graph g2 = new Graph(graphX,graphY,simWidth,simWidth,sim1.outer.values,c2,0,maxTemp);
    Graph g3 = new Graph(graphX,graphY,simWidth,simWidth,sim1.differences,c3,0,maxTemp);
    
    //sim2
    graphX+=simWidth;
    Graph g4 = new Graph(graphX,graphY,simWidth,simWidth,sim2.inner.values,c1,0,maxTemp);
    Graph g5 = new Graph(graphX,graphY,simWidth,simWidth,sim2.outer.values,c2,0,maxTemp);
    Graph g6 = new Graph(graphX,graphY,simWidth,simWidth,sim2.differences,c3,0,maxTemp);
    
    //sim3
    graphX+=simWidth;
    Graph g7 = new Graph(graphX,graphY,simWidth,simWidth,sim3.inner.values,c1,0,maxTemp);
    Graph g8 = new Graph(graphX,graphY,simWidth,simWidth,sim3.outer.values,c2,0,maxTemp);
    Graph g9 = new Graph(graphX,graphY,simWidth,simWidth,sim3.differences,c3,0,maxTemp);
    
    //sim4
    graphX+=simWidth;
    Graph g10 = new Graph(graphX,graphY,simWidth,simWidth,sim4.inner.values,c1,0,maxTemp);
    Graph g11 = new Graph(graphX,graphY,simWidth,simWidth,sim4.outer.values,c2,0,maxTemp);
    Graph g12 = new Graph(graphX,graphY,simWidth,simWidth,sim4.differences,c3,0,maxTemp);
    
    //draw graphs to screen with dividers in between
    g1.draw();
    g2.draw();
    g3.draw();
    
    stroke(divider);
    line(simWidth,yOffset,simWidth,height);
    
    g4.draw();
    g5.draw();
    g6.draw();
    
    stroke(divider);
    line(simWidth*2,yOffset,simWidth*2,height);
    
    g7.draw();
    g8.draw();
    g9.draw();
    
    stroke(divider);
    line(simWidth*3,yOffset,simWidth*3,height);
    
    g10.draw();
    g11.draw();
    g12.draw();
    
    noStroke();
    
    //draw legend onto screen using above specified colours
    fill(c1);
    rect( 5, simWidth+yOffset-36,6,11);
    fill(0);
    text("Inner temperature",15,simWidth+yOffset-24);
    
    fill(c2);
    rect( 5,simWidth+yOffset-24,6,11);
    fill(0);
    text("Outer temperature",15,simWidth+yOffset-12);
    
    fill(c3);
    rect( 5,simWidth+yOffset-12,6,11);
    fill(0);
    text("Outer-Inner temperature difference",15,simWidth+yOffset);
    
    fill(100);
    
    //draw heat gradient onto screen
    image(heatGrad,width-simWidth,simWidth+yOffset-10);
    
    //draw textual information above simulations
    text(minTemp+" C", width-simWidth,simWidth+yOffset-11);
    text(maxTemp+" C", width-simWidth+255-50,simWidth+yOffset-11);
    text("Object Temperature Scale",width-simWidth+40,simWidth+yOffset-24);
    
    
    text(minTemp+" C",0,height-5);
    text(maxTemp+" C",0,simWidth+12+yOffset);
    
    //Draw title bar
    text(time,5,10);
    fill(0);
    text("Heat Dispersion in 3D Shapes in an Oven",width/2-100,15);
  }
}


void keyPressed() {
  //pause simulation on space bar
  if (key == ' ') {
    paused = !paused;
  }
}


/*
 The class inherit all the fields, constructors and functions 
 of the java.awt.Polygon class, including contains(), xpoint,ypoint,npoint
 
 This is simply to make drawing polygons much more practical in our simulation
 */
class Poly extends java.awt.Polygon {

    public Poly(int[] x, int[] y, int n) {
        //call the java.awt.Polygon constructor
        super(x, y, n);
    }

    void drawMe() {
        beginShape();
        for (int i = 0; i < npoints; i++) {
            vertex(xpoints[i], ypoints[i]);
        }
        endShape(CLOSE);
    }
}


/*
This class creates and draws graphs of time-series data
*/

class Graph{

  //local locatoin and size variables  
  int x;
  int y;
  int w;
  int h;
  float mx;
  float mn;
  
  //colour of this graph
  color f;
  
  //graph data
  ArrayList<Float> series;
  
  public Graph(int x, int y, int w, int h, ArrayList<Float> data,color f,float mn,float mx){
    this.x =x;
    this.w=w;
    this.y=y;
    this.h=h;
    this.mn=mn;
    this.mx=mx;
    series=data;
    this.f = f;
  }
  
  public void draw(){
       
    //calculate differences between drawing points
    float difference = series.size()*(1.0/w);
    //calculate height scale
    float scaleFactor = h*1.0/(mx-mn);
    
    //init drawing properties
    strokeWeight(1);
    stroke(f);
    fill(f);
    
    //finally iterate through pixels and draw data
    for(int i=0;i<w;i++){
      
      if( i*difference < series.size()-1 && (i+1)*difference < series.size()-1){
        line( x+i,y+h-series.get((int)((i)*difference))*scaleFactor,
              x+i+1,y+h-series.get((int)((i+1)*difference))*scaleFactor
        );
      }
    }
  }
  
}


/*

Thermometer class, which holds time-series data

In the future this should use temporary files to better manage memory,
but under constraints of competition this was never necessary on my 
machine, and so I didn't use the time to implement the improvements.

*/
class Thermometer {
    
    //data
    ArrayList<Float> values;
    
    //thermometer location in mesh
    int x;
    int y;
    int z;

    public Thermometer(int x, int y, int z) {
        this.x = x;
        this.y = y;
        this.z = z;
        values = new ArrayList<Float>();
    }

    //update this thermometer from the mesh
    public void update(float[][][] nodes) {
      try{
        float val = nodes[x][y][z];
        values.add((Float)val);
      }
      catch(Exception e){
        print(e);
      }
    }
    
    //calculate area under curve, i.e. total heat absorbed
    public float calculateTotalHeat(){
      
      float heat = 0;
      for(int i=0;i<values.size();i++){
        heat += values.get(i);
      }
      
      return heat;
      
    }
}


/*
This class runs an individual simulation.

*/
class Simulation {

    //additional data
    float currentTime = 0;
  
    //polygon info
    int polygonSides = 4;
    float polyRadius = 50;
    int cx=w/2;
    int cy=l/2;
  
    int cornerX = cx;
    int cornerY = cy;

    int centerX = cx;
    int centerY = cy;
    
    //simulation coordinates on screen
    int drawX;
    int drawY;
    
    //diagnostic info
    float maxDifferenceRatio = 0;
    
    //required variables
    float [][][] nodes = new float[w][l][h];
    ArrayList<int[]> object = new ArrayList();
    ArrayList<int[]> objectEdges = new ArrayList();
    ArrayList<Float> differences = new ArrayList();
    Polygon surface = new Polygon();
    
    //thermometers
    Thermometer inner,outer;  
    
    //initiallize from polygon
    public Simulation(Polygon poly, int x, int y){
        
        surface = poly;

        drawX = x;
        drawY = y;

        centerX = center(surface)[0];
        centerY = center(surface)[1];
        cornerX = surface.xpoints[0];
        cornerY = surface.ypoints[0];

        inner = new Thermometer(centerX,centerY,h/2);
        outer = new Thermometer(cornerX,cornerY+1,h/2);

        //generate heat distribution and mesh
        for (int i = 0; i < nodes.length; i++) {
            for (int j = 0; j < nodes[0].length; j++) {

                for (int k = 0; k < nodes[0][0].length; k++) {
                    nodes[i][j][k] = airTemp;

                    if (surface.contains(i, j)) {
                        nodes[i][j][k] = surfaceTemp;
                        int[] loc = {
                            i, j, k
                        };
                        object.add(loc);

                        if (edge(i, j, k, nodes, surface)) {
                            objectEdges.add(loc);
                        }
                    }
                }
            }
        }
        
        
        
        PImage img = mapToImage(nodes);
        image(img, 0, 0);
        text(currentTime, 10, 20);
      
    }
    
    //initialize (regular polygon) from side count and length
    //TODO: Refactor this to a cleaner form
    public Simulation(int polygonSides,int x,int y,float polyLen) {
        
        this.polygonSides = polygonSides;
      
        drawX = x;
        drawY = y;

        initPolygon(polyLen);
        
        inner = new Thermometer(centerX,centerY,h/2);
        outer = new Thermometer(cornerX+2,cornerY+2,h/2);

        if(polygonSides == 3){
          outer.x+=2;
        }
  
        //generate heat distribution
        for (int i = 0; i < nodes.length; i++) {
            for (int j = 0; j < nodes[0].length; j++) {

                for (int k = 0; k < nodes[0][0].length; k++) {
                    nodes[i][j][k] = airTemp;

                    if (surface.contains(i, j)) {
                        nodes[i][j][k] = surfaceTemp;
                        int[] loc = {
                            i, j, k
                        };
                        object.add(loc);

                        if (edge(i, j, k, nodes, surface)) {
                            objectEdges.add(loc);
                        }
                    }
                }
            }
        }

        PImage img = mapToImage(nodes);
        image(img, 0, 0);
        text(currentTime, 10, 20);
    }
    
    //create a polygon meeting initialized specifications
    public void initPolygon(float polyLen){
      //initialize polygon
        
        cx-=polyRadius/2;
        cy-=polyRadius/2;
        
        float px = cx;
        float py = cy;
        
        polyRadius = polyLen;
        
        for (int i = 0; i < polygonSides; i++) {
            //Note, I used the Java Math sin and cosine here
            //it doesn't make a difference
            px += polyRadius * Math.cos(2*PI * 1.0 / polygonSides * i);
            py += polyRadius * Math.sin(2*PI * 1.0 / polygonSides * i);
            surface.addPoint((int)px, (int)py);
        }
        
        centerX = center(surface)[0];
        centerY = center(surface)[1];
        cornerX = cx;
        cornerY = cy;
    }
    
    //check if a mesh voxel is on an edge
    boolean edge(int x, int y, int z, float[][][] mesh, Polygon surface) {

        return ((surface.contains(x - 1, y)
                && surface.contains(x + 1, y)
                && surface.contains(x, y - 1)
                && surface.contains(x, y + 1)
                && z > 0 //&& z < mesh[x][y].length-1
                ) == false
                && surface.contains(x, y));
    }
    
    //draw simulation
    void draw() {
        //initialize drawing context
        PGraphics graph = createGraphics(simWidth,simWidth);

        //apply heating and calculations
        nodes = addSurfaceHeat(nodes);
        nodes = updateTemperature(nodes);
        
        //update thermometers
        inner.update(nodes);
        outer.update(nodes);
        differences.add( abs(inner.values.get(inner.values.size()-1) - outer.values.get(outer.values.size()-1)));
        
        //draw image
        PImage img = mapToImage(nodes);

        graph.beginDraw();

        graph.image(img, 0, 0);
        graph.fill(0);
        //graph.text(currentTime, 10, 20);

        graph.stroke(0);
        graph.noFill();
        
        
        //draw surface of object
        graph.beginShape();
        for (int i = 0; i < surface.npoints; i++) {
            graph.vertex(surface.xpoints[i], surface.ypoints[i]);
        }
        graph.endShape(CLOSE);
        
        
        graph.fill(0, 100, 255);
        graph.stroke(255);
        graph.strokeWeight(0.5);
        
        //draw thermometer points
        graph.ellipse(centerX, centerY, 4, 4);
        graph.ellipse(cornerX, cornerY, 4, 4);
        
        fill(100);
        
        //write text for simulation information
        graph.text("Center temperature:\t" + inner.values.get(inner.values.size()-1), 10, 17);
        graph.text("Edge temperature:\t" + outer.values.get(outer.values.size()-1), 10, 34);
        graph.text("Air temperature:\t" + airTemp, 10, 51);
        graph.text("Mesh Depth:\t" + h,10,68);
        
        //float totalHeatDifference = 0;
        //for(int i=0;i<differences.size();i++) totalHeatDifference+= differences.get(i);
        //graph.text("Total Heat delta:\t"+ totalHeatDifference,10,108);
        
        float heatingDelta = outer.calculateTotalHeat()/inner.calculateTotalHeat();
        maxDifferenceRatio = max(maxDifferenceRatio,heatingDelta);
        graph.text("Outer to Inner heating diff. ratio:\t"+heatingDelta,10,85 );
        graph.text("Max heating difference ratio:\t"+maxDifferenceRatio,10,102);
        graph.endDraw();

        //finally draw to screen
        image(graph,drawX,drawY);
        //  exit();
    }
    
    //find center of polygon
    int[] center(Polygon p) {
        int x = 0;
        int y = 0;

        for (int i = 0; i < p.npoints; i++) {
            x += p.xpoints[i];
            y += p.ypoints[i];
        }

        x /= p.npoints;
        y /= p.npoints;

        int[] loc = {
            x, y
        };

        return loc;
    }
    
    //get voxel temperature
    float getNodeTemperature(int x, int y, int z) {
        return objectContains(x, y, z) ? nodes[x][y][z] : -1;
    }
    
    //check if voxel mesh contains a point
    //note this only works because we are using basic polygonal prisms
    boolean objectContains(int x, int y, int z) {
        return surface.contains(x, y);
    }

    //update temperature within surface
    float[][][] updateTemperature(float[][][] nodeList) {

        float[][][] newNodeList = nodeList;
       
        for (int i = 0; i < object.size(); i++) {
            int x = object.get(i)[0];
            int y = object.get(i)[1];
            int z = object.get(i)[2];
            //apply heat diffusion
            newNodeList[x][y][z] = calcTemp(nodeList, x, y, z);
            //    print(x+" "+y+"\t");
        }
        
        return newNodeList;
    }
    
    //This method was skipped in practice since under time constraint
    //it was sufficient to simulate for aluminum
    float calculateInternalEnergy(float row, float cp, float T1, float T2) {
        return row * cp * (T2 - T1);
    }
    //same as above
    float calculateDiffusivity() {
        //return 84.14;
        //return 8.418E-5;
        return k / row / C;
    }
    
    //same as above
    float meshFourierNumber(float diffusivity, float deltaT, float deltaX) {
        return 1;
        //return diffusivity*deltaT/pow(deltaX, 2);
    }
    
    //same as above
    float calculateHeatDifference(float h) {
        return 0;
    }
    
    //Apply heat from surrounding convection oven to edges
    //note that the top of the pan was assumed to absorb heat from the
    //brownie, hence it didn't need to receive convection heating
    //we made this assumption because metal heats faster than
    //brownie mix
    float[][][] addSurfaceHeat(float[][][] list) {

        for (int i = 0; i < objectEdges.size(); i++) {
            int x = objectEdges.get(i)[0];
            int y = objectEdges.get(i)[1];
            int z = objectEdges.get(i)[2];
            list[x][y][z] += (airTemp - list[x][y][z]) * delx * dely * heatTransfer;
        }
        return list;
    }
    
    //check if point is inside surface
    boolean inSurface(int x, int y, int z) {
        int[] loc = {
            x, y, z
        };
        return object.contains(loc);
    }
    
    //better method of surface checking; works for z
    boolean inSurface2(int x, int y, int z, float list[][][], Polygon surface) {
        return surface.contains(x, y) && z >= 0 && z < list[0][0].length;
    }

    //calculate desired temperature
    float calcTemp(float[][][] list, int x, int y, int z) {
        //calculated using an explicit method
        //new temperature node
        
        //applying heat equation
        float Tn = list[x][y][z];

        float Tz = 0;
        float Ty = 0;
        float Tx = 0;

        Tx += inSurface2(x + 1, y, z, list, surface) ? list[x + 1][y][z] : list[x][y][z];
        Tx += inSurface2(x - 1, y, z, list, surface) ? list[x - 1][y][z] : list[x][y][z];

        Ty += inSurface2(x, y + 1, z, list, surface) ? list[x][y + 1][z] : list[x][y][z];
        Ty += inSurface2(x, y - 1, z, list, surface) ? list[x][y - 1][z] : list[x][y][z];

        Tz += inSurface2(x, y, z - 1, list, surface) ? list[x][y][z - 1] : list[x][y][z];
        Tz += inSurface2(x, y, z + 1, list, surface) ? list[x][y][z + 1] : list[x][y][z];

        Tn = (Tx + Ty + Tz) / 6;
        return Tn;

    }

    //converts an integer map to a PImage
    PImage mapToImage(float[][][] img) {
        //generate final image
        PImage finalImage = new PImage(img.length, img[0].length);

        //paste data into final image
        for (int x = 0; x < img.length; x++) {
            for (int y = 0; y < img[x].length; y++) {
                float val = 0;

                for (int z = 0; z < img[x][y].length; z++) {
                    val += img[x][y][z];
                }

                float heat = 255 / (maxTemp - minTemp) * val / img[0][0].length;
                if (val < 0) {
                    finalImage.set(x, y, color(0));
                } else {
                    finalImage.set(x, y, color(heat));
                }
            }
        }

        currentTime += timestep;
        //return the final image map
        return finalImage;
    }
}

