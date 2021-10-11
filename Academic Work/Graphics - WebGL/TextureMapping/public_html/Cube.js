"use strict";

function getTextureCoords (horizontalSegments, verticalSegments) {
  let textureCoords = [];
  for (let i = 0; i < verticalSegments; i++) {
    for (let j = 0; j < horizontalSegments; j++) {
      textureCoords = [
        ...textureCoords,
        0,0,1,0,1,1,0,1
//        j / horizontalSegments, i / verticalSegments,
//        (j + 1) / horizontalSegments, i / verticalSegments,
//        (j + 1) / horizontalSegments, (i + 1) / verticalSegments,
//        j / horizontalSegments, (i + 1) / verticalSegments
      ];
    }
    return textureCoords;
  }}
  
function getIndices (horizontalSegments, verticalSegments) {
  let indeces = [];
  for (let i = 0; i < horizontalSegments * verticalSegments; i += 4) {
    indeces = [
      ...indeces,
      i, i + 1, i + 2,
      i, i + 2, i + 3
    ];
  }
}

function getNormals (vertices) {
  let normals = [];
  for (let i = 0; i < vertices/4; i += 4) {
    let normalX = vertices[i] + vertices[i + 3] + vertices[i + 6] + vertices[i + 9];
    let normalY = vertices[i + 1] + vertices[i + 4] + vertices[i + 7] + vertices[i + 10];
    let normalZ = vertices[i + 2] + vertices[i + 5] + vertices[i + 8] + vertices[i + 11];
    const hypot = Math.hypot(normalX, normalY, normalZ);

    normalX /= hypot;
    normalY /= hypot;
    normalZ /= hypot;

    for (let i = 0; i < 4; i++) {
      normals = [
        ...normals,
        normalX,
        normalY,
        normalZ
      ];
    }
  }
}
  
function getPointFromAngles(horizontal, vertical, scale=1) {
  return [
    Math.cos(vertical) * Math.cos(horizontal) * scale,
    Math.sin(vertical) * scale,
    Math.cos(vertical) * Math.sin(horizontal) * scale
  ];
}

function getUVSphere (hDivs, vDivs, scale=1) {
  var vertices = [];
  var hInterval = 2 * Math.PI / hDivs;
  var vInterval = Math.PI / vDivs;
  var hAngle = 0;
  var vAngle = Math.PI / 2;
  for (let i = 0; i < vDivs; i++) {
    for (let j = 0; j < hDivs; j++) {
      vertices.push(...getPointFromAngles(hAngle, vAngle + vInterval, scale));
      vertices.push(...getPointFromAngles(hAngle + hInterval, vAngle + vInterval, scale));
      vertices.push(...getPointFromAngles(hAngle, vAngle, scale));
       vertices.push(...getPointFromAngles(hAngle + hInterval, vAngle, scale));
      
      hAngle += hInterval;
    }
    vAngle -= vInterval;
  }
  const normals = getNormals(vertices);
  const textCoords = getTextureCoords(hDivs, vDivs);
  const indices = getIndices(hDivs, vDivs);
  
  //return vertices;
  return {vertices, normals, textCoords, indices};
}

// The Cube class makes a unit cube at the center of the world
function Cube(gl, shaderProgram) {
    var t = this;
    this.gl = gl;

    // Create a Vertex Array Object.  This remembers the buffer bindings and
    // vertexAttribPointer settings so that we can reinstate them all just using
    // bindVertexArray.
    this.vao = gl.createVertexArray();  // create and get identifier
    gl.bindVertexArray(this.vao);

    this.numVertices = 400;  // total number of vertices
    this.faceVertices = 4;  // number of vertices in a face

    const UVSphere = getUVSphere(10,10,1);
    console.log(UVSphere);
    const verts = UVSphere.vertices;
    const texCoords = UVSphere.textCoords;
//    var verts = [0.5, -0.5, -0.5, // +x face
//        0.5, 0.5, -0.5,
//        0.5, -0.5, 0.5,
//        0.5, 0.5, 0.5,
//
//        -0.5, -0.5, 0.5, // -x face
//        -0.5, 0.5, 0.5,
//        -0.5, -0.5, -0.5,
//        -0.5, 0.5, -0.5,
//
//        -0.5, 0.5, -0.5, // +y face
//        -0.5, 0.5, 0.5,
//        0.5, 0.5, -0.5,
//        0.5, 0.5, 0.5,
//
//        0.5, -0.5, -0.5, // -y face
//        0.5, -0.5, 0.5,
//        -0.5, -0.5, -0.5,
//        -0.5, -0.5, 0.5,
//
//        -0.5, -0.5, 0.5, // +z face
//         0.5, -0.5, 0.5,
//        -0.5, 0.5, 0.5,
//        0.5, 0.5, 0.5,
//
//        -0.5, 0.5, -0.5, // -z face
//         0.5, 0.5, -0.5,
//        -0.5, -0.5, -0.5,
//        0.5, -0.5, -0.5
//    ];

var colors = [];
for(var i = 0; i < 100; i++)
{
    colors = [...colors,1,1,1,1,1,1,1,1,1,1,1,1];
}
        
            //[1, 0, 0,   // red
//        1, 0, 0,
//        1, 0, 0,
//        1, 0, 0,
//
//        0, 1, 0,    // green
//        0, 1, 0,
//        0, 1, 0,
//        0, 1, 0,
//
//        0, 0, 1,     // blue
//        0, 0, 1,
//        0, 0, 1,
//        0, 0, 1,
//
//        1, 1, 0,      // yellow
//        1, 1, 0,
//        1, 1, 0,
//        1, 1, 0,
//
//        1, 0, 1,      // magenta
//        1, 0, 1,
//        1, 0, 1,
//        1, 0, 1,
//
//        0, 1, 1,      // cyan
//        0, 1, 1,
//        0, 1, 1,
//        0, 1, 1
//    ];

//
//    var texCoords = 
//       [0.5, 0,           // + x
//        0.5, 0.5,
//        0, 0,
//        0, 0.5,
//        
//        0.5, 1,            // -x  -> 90 degree 
//        0, 0.5,
//        1, 0.5,
//        0.5, 1,
//
//        1, 0,           // + y
//        1, 1,
//        0, 0,
//        0, 1,
//
//        0, 0,           // - y
//        0, 1,
//        1, 0,
//        1, 1,
//
//        0, 2,            // +z   -> small
//        2, 2,            
//        0, 0,           
//        2, 0,
//
//        0, 1,            // -z   
//        1, 1,            
//        0, 0,           
//        1, 0
//    ];



    var floatBytes = 4;  // number of bytes in a float value

    // Create and load the vertex positions
    this.cubeVertVB = gl.createBuffer();  // get unique buffer ID number
    gl.bindBuffer(gl.ARRAY_BUFFER, this.cubeVertVB);
    gl.bufferData(gl.ARRAY_BUFFER, flatten(verts), gl.STATIC_DRAW);
    this.vPosition = gl.getAttribLocation(shaderProgram, "vPosition");
    gl.vertexAttribPointer(this.vPosition, 3, gl.FLOAT, false, 3 * floatBytes, 0);
    gl.enableVertexAttribArray(this.vPosition);


    // Create and load the vertex colors
    this.cubeColorVB = gl.createBuffer();  // get unique buffer ID number
    gl.bindBuffer(gl.ARRAY_BUFFER, this.cubeColorVB);
    gl.bufferData(gl.ARRAY_BUFFER, flatten(colors), gl.STATIC_DRAW);
    this.vColor = gl.getAttribLocation(shaderProgram, "vColor");
    gl.vertexAttribPointer(this.vColor, 3, gl.FLOAT, false, 3 * floatBytes, 0);
    gl.enableVertexAttribArray(this.vColor);

    // Create and load the vertex texture coords
    this.cubeTextureVB = gl.createBuffer();  // get unique buffer ID number
    gl.bindBuffer(gl.ARRAY_BUFFER, this.cubeTextureVB);
    gl.bufferData(gl.ARRAY_BUFFER, flatten(texCoords), gl.STATIC_DRAW);
    this.vTexCoord = gl.getAttribLocation(shaderProgram, "vTexCoord");
    gl.vertexAttribPointer(this.vTexCoord, 2, gl.FLOAT, false, 2 * floatBytes, 0);
    gl.enableVertexAttribArray(this.vTexCoord);


    // Get uniform variable location for transform matrices
    this.projectionMat = gl.getUniformLocation(shaderProgram, "projectionMat");
    this.viewMat = gl.getUniformLocation(shaderProgram, "viewMat");
    this.modelMat = gl.getUniformLocation(shaderProgram, "modelMat");


    // Texture uniform in frag shader
    this.fTexSampler = gl.getUniformLocation(shaderProgram, "fTexSampler");
    this.fShowTexture = gl.getUniformLocation(shaderProgram, "fShowTexture");
    this.fShowColor = gl.getUniformLocation(shaderProgram, "fShowColor");

    // Init the textureMap - this starts the loading and once loaded sets
    //   some initial paramters.  It also sets a variable to let us know
    //   when it is done loading since it loads asynch.
    this.InitTexture("brick.png");
    

    gl.bindVertexArray(null);  // un-bind our vao

}

// Render function that draws the cube
// We need the 3 transformations matrices
//   as well as info about colors, textures and texture filters
Cube.prototype.Render = function (projectionMat, viewMat, modelMat, colorsOn, textureOn) {
    var gl = this.gl;

    // Set up buffer bindings and vertex attributes
    // We did this with the VAO so just use it
    gl.bindVertexArray(this.vao);

    // Set transformation matrices for shader
    gl.uniformMatrix4fv(this.projectionMat, false, flatten(projectionMat));
    gl.uniformMatrix4fv(this.viewMat, false, flatten(viewMat));
    gl.uniformMatrix4fv(this.modelMat, false, flatten(modelMat));

  
    // Set up texture
    // Note that for this lab we will always have a texture - either the image
    //   or white.  Given the booleans passed to the shader, the white texture
    //   won't even be applied if no texture is selected.  But the small white
    //   texture prevents shader warnings that no active texture is bound.
    //   We could have made two differrent shaders to do this as well
    
    gl.activeTexture(gl.TEXTURE0);  // which of the multiple texture units to use
    gl.uniform1i(this.fTexSampler, 0); // The texture unit to use
    
    if (textureOn && this.textureLoaded) {
       
        gl.bindTexture(gl.TEXTURE_2D, this.texture);  // The image
        gl.uniform1i(this.fShowColor, colorsOn);
        gl.uniform1i(this.fShowTexture, 1);
        gl.texParameteri(gl.TEXTURE_2D, gl.TEXTURE_MAG_FILTER, gl.NEAREST);
        gl.texParameteri(gl.TEXTURE_2D, gl.TEXTURE_MIN_FILTER, gl.NEAREST);
        
    } else {
        gl.bindTexture(gl.TEXTURE_2D, this.whiteTexture);  // The white texture
        gl.uniform1i(this.fShowColor, colorsOn);
        gl.uniform1i(this.fShowTexture, 0);
    }

    
    // Draw each face as a 2-triangle (4-vertex) strip
    for (var start = 0; start < this.numVertices; start += this.faceVertices) {
        gl.drawArrays(gl.TRIANGLE_STRIP, start, this.faceVertices);
    }
    
    gl.bindVertexArray(null);  // un-bind our vao

};


Cube.prototype.InitTexture = function (textureURL) {
    var gl = this.gl;
    
    // First make a white texture for when we don't want to have a texture
    //   This prevents shader warnings even if we don't sample from it
    this.whiteTexture = gl.createTexture();
    gl.bindTexture(gl.TEXTURE_2D, this.whiteTexture);
    gl.texImage2D(gl.TEXTURE_2D, 0, gl.RGBA, 1, 1, 0,
                  gl.RGBA, gl.UNSIGNED_BYTE, new Uint8Array([255, 255, 255, 255]));

    // Load the texture from url (with generated mipmaps)
    this.textureLoaded = false;

    var texture = this.texture = gl.createTexture();
    var textureImage = new Image();
    var t = this;

    // Set up function to run asynchronously after texture image loads
    textureImage.onload = function () {
        gl.activeTexture(gl.TEXTURE0);
        gl.bindTexture(gl.TEXTURE_2D, texture);
        gl.pixelStorei(gl.UNPACK_FLIP_Y_WEBGL, true);
        gl.texImage2D(gl.TEXTURE_2D, 0, gl.RGBA, gl.RGBA, gl.UNSIGNED_BYTE, textureImage);
      
        gl.texParameteri(gl.TEXTURE_2D, gl.TEXTURE_WRAP_S, gl.REPEAT);
        gl.texParameteri(gl.TEXTURE_2D, gl.TEXTURE_WRAP_T, gl.REPEAT);
        
        //Determine min and mag based on radio buttons
        //Mag Nearestrameteri(gl.TEXTURE_2D, gl.TEXTURE_MAG_FILTER, gl.NEAREST
        if(document.getElementById("magNearest").checked){gl.texParameteri(gl.TEXTURE_2D, gl.TEXTURE_MAG_FILTER, gl.NEAREST);}
        //Mag Linear
        else if(document.getElementById("magLinear").checked){gl.texParameteri(gl.TEXTURE_2D, gl.TEXTURE_MAG_FILTER, gl.LINEAR);}
        //Min Nearest
        if(document.getElementById("minNearest").checked){gl.texParameteri(gl.TEXTURE_2D, gl.TEXTURE_MIN_FILTER, gl.NEAREST);}
        //Min Linear
        else if(document.getElementById("minLinear").checked){ gl.texParameteri(gl.TEXTURE_2D, gl.TEXTURE_MIN_FILTER, gl.LINEAR);}
        //Min LML
        else if(document.getElementById("minLinearMipLinear").checked){gl.texParameteri(gl.TEXTURE_2D, gl.TEXTURE_MIN_FILTER, gl.LINEAR_MIPMAP_LINEAR);}
      
        gl.generateMipmap(gl.TEXTURE_2D);  // incase we need min mipmap
        
        t.textureLoaded = true;  // flag texture load complete
    };

    textureImage.src = textureURL;  // start load of texture image
    
};