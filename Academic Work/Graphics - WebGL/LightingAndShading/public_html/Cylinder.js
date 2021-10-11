"use strict";

class Foo {
      constructor(horizontalDiv, verticalDiv) {
        this.horizontalDiv = horizontalDiv;
        this.verticalDiv = verticalDiv;

        this.vertices = [];
        this.indices = [];
        this.normals = [];
        this.texCoords = [];

        let horizontalAng = 0;
        let horizontalDiff = (2 * Math.PI) / this.horizontalDiv;

        let verticalAng = Math.PI / 2;
        let verticalDiff = -(Math.PI / this.verticalDiv);


        for (var i = 0; i <= this.verticalDiv; i++) {

          let cosVert = Math.cos(verticalAng);
          let sinVert = Math.sin(verticalAng);

          for (var j = 0; j <= this.horizontalDiv; j++) {

            let cosHor = Math.cos(horizontalAng);
            let sinHor = Math.sin(horizontalAng);

            // z = (r * cos(verticalAng)) * cos(horizontalAng)
            // x = (r * cos(verticalAng)) * sin(horizontalAng)
            // y = r * sin(veritcalAng)

            let x = cosVert * sinHor;
            let y = sinVert;
            let z = cosVert * cosHor;

            this.vertices.push(x, y, z);
            this.normals.push(x, y, z);

            this.texCoords.push(j / this.horizontalDiv);
            this.texCoords.push(i / this.verticalDiv);

            horizontalAng += horizontalDiff;
          }

          verticalAng += verticalDiff;
        }


        for (var i = 0; i < this.verticalDiv; i++) {
          let k1 = i * (this.horizontalDiv + 1);
          let k2 = k1 + this.horizontalDiv + 1;

          for (var j = 0; j < this.horizontalDiv; j++) {

            if (i != 0) {
              this.indices.push(k1);
              this.indices.push(k2);
              this.indices.push(k1 + 1);
            }

            if (i != (this.verticalDiv - 1)) {
              this.indices.push(k1 + 1);
              this.indices.push(k2);
              this.indices.push(k2 + 1);
            }

            k1++;
            k2++;
          }
        }
      }
    }

// The Cylinder class that defines and draws an open cylinder
// We take in the scale in the constructor to modify the vertex buffer points
//   by the scale on construction - so we never need to scale each point when rendering
function Cylinder(gl, shaderProgram, scale) {
    var t = this;
    this.gl = gl;
    
    // Create a Vertex Array Object.  This remembers the buffer bindings and
    // vertexAttribPointer settings so that we can reinstate them all just using
    // bindVertexArray.
    this.vao = gl.createVertexArray();  // create and get identifier
    gl.bindVertexArray(this.vao);


    // A cylinder is just a triangle strip alternating top and bottom circle points
    var radius = 1.0 * scale;
    var inRad = 0.5 * scale;
    this.numSections = 100;
   
    // Create separate CPU arrays for vertex points and normals
    // Note we do not create colors as they will be computed in the shaderss
    var verts = [];
    var normals = [];
//    for (var i=0; i<=this.numSections; i++) {  // we go <= to get the circle to close
//       var angleDeg = (360.0 / this.numSections)*i;
//       var angleRad = angleDeg*Math.PI/180.0;
//       var x = radius * Math.cos(angleRad);
//       var z = radius * Math.sin(angleRad);
//       var y = 1.0 * scale;
//       verts.push(vec3(x, y, z));   // top circle point
//       normals.push(normalize(vec3(x, 0.0, z)));  // unit vector
//       
//       verts.push(vec3(x, -y, z));  // bottom circle point
//       normals.push(normalize(vec3(x, 0.0, z)));  // unit vector
//       
//    }
    
    const f = new Foo(10,10);
    verts = f.vertices;
    normals = f.normals;
    var floatBytes = 4;  // number of bytes in a float value
        
    // Create and load the vertex positions
    this.cylVertVB = gl.createBuffer();  // get unique buffer ID number
    gl.bindBuffer(gl.ARRAY_BUFFER, this.cylVertVB );
    gl.bufferData(gl.ARRAY_BUFFER, flatten(verts), gl.STATIC_DRAW );
    this.vPosition = gl.getAttribLocation(shaderProgram, "vPosition");
    gl.vertexAttribPointer(this.vPosition, 3, gl.FLOAT, false, 3 * floatBytes, 0);
    gl.enableVertexAttribArray(this.vPosition);
     
    
    // Create and load the vertex normals
    this.cylNormalVB = gl.createBuffer();  // get unique buffer ID number
    gl.bindBuffer(gl.ARRAY_BUFFER, this.cylNormalVB );
    gl.bufferData(gl.ARRAY_BUFFER, flatten(normals), gl.STATIC_DRAW );
    this.vNormal = gl.getAttribLocation(shaderProgram, "vNormal");
    gl.vertexAttribPointer(this.vNormal, 3, gl.FLOAT, false, 3 * floatBytes, 0);
    gl.enableVertexAttribArray(this.vNormal);
    
    
    // Get uniform variable location for transform matrices
    // Reall they need to be separately due to lighting calculations
    this.projectionMat = gl.getUniformLocation(shaderProgram, "projectionMat");
    this.viewMat = gl.getUniformLocation(shaderProgram, "viewMat");
    this.modelMat = gl.getUniformLocation(shaderProgram, "modelMat");
    
    // Get uniform variable locations for light info
    // At the moment, we only send one color and use it for ambient, diffuse and specular
    this.lightPosition = gl.getUniformLocation(shaderProgram, "lightPosition");
    this.lightColor = gl.getUniformLocation(shaderProgram, "lightColor");
    this.ambientFactor = gl.getUniformLocation(shaderProgram, "ambientFactor");
    
    // Get uniform variable locations for material properties (K)
    // At the moment, we only send one color and use it for ambient, diffuse and specular
    this.materialColor = gl.getUniformLocation(shaderProgram, "materialColor");
    this.materialShiny = gl.getUniformLocation(shaderProgram, "shiny");
    
    // Get uniform variable locations for Phong components to use
    // Basically these are booleans to control which components the shader uses
    // Note that we could have just created 3 different shaders and then
    //  selected the correct shader based on these boolens.  However, I feel this
    //  way is less code and a bit easier to understand.
    this.ambientOn = gl.getUniformLocation(shaderProgram, "ambientOn");
    this.diffuseOn = gl.getUniformLocation(shaderProgram, "diffuseOn");
    this.specularOn = gl.getUniformLocation(shaderProgram, "specularOn");
   
    
    gl.bindVertexArray(null);  // un-bind our vao
   
}

// Render function that draws the cylinder
// We need the 3 transformations matrices as well as the light position
// We are also passed 4 booleans to control the redering details
Cylinder.prototype.Render = function(projectionMat, viewMat, modelMat, lightPosition, showEdges, ambientOn, diffuseOn, specularOn) {
   var gl = this.gl; 
    
   // Set up buffer bindings and vertex attributes
   // We did this with the VAO so just use it
   gl.bindVertexArray(this.vao);
  
   // Set transformation matrices for shader
   gl.uniformMatrix4fv(this.projectionMat, false, flatten(projectionMat));
   gl.uniformMatrix4fv(this.viewMat, false, flatten(viewMat));
   gl.uniformMatrix4fv(this.modelMat, false, flatten(modelMat));
   
  
   // Light color - hard-coded here to be white
   // Note currently this color is used for ambent, diffuse and specular
   var lColor = vec3(1.0, 1.0, 1.0);
   var ambientFactor = 0.3;  // 30% ambient on everything
       
   // Pass in the light info
   // No false for these as 2nd parameter
   gl.uniform3fv(this.lightPosition, flatten(lightPosition));
   gl.uniform3fv(this.lightColor, flatten(lColor));
   gl.uniform1f(this.ambientFactor, ambientFactor);
    
   
   // Pass in the material color - hard-coded here to be semi-shiny blue
   // Note currently this color is used for ambent, diffuse and specular
   var mColor = vec3(0.8, 0.1, 0.1);
   var mShiny = 50.0;
   gl.uniform3fv(this.materialColor, flatten(mColor));
   gl.uniform1f(this.materialShiny, mShiny);
   
   
   // Pass in which components to use (booleans pass as 0/1 ints)
   gl.uniform1i(this.ambientOn, ambientOn);
   gl.uniform1i(this.diffuseOn, diffuseOn);
   gl.uniform1i(this.specularOn, specularOn);

   
   // Tell the pipeline to draw the triangles or wireframe
   if (!showEdges) {   
       for(var i = 0; i < 10; i++)
       {
           gl.drawArrays(gl.TRIANGLE_STRIP, i*this.numSections, this.numSections*2+2);
       }
   } else {
       for(var i = 0; i < this.numSections * 2; i += 1)
       {
            gl.drawArrays(gl.LINE_LOOP, i, 12);
       }
   }
 
          
  gl.bindVertexArray(null);  // un-bind our vao
 
};