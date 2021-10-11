/**
 * Lab 2 - COMP3801 Spring 2021
 *   Basic WebGL2 shaders, mouse events and coordinates
 */

"use strict";

// Constructor
//
// @param canvasID - string containing name of canvas to render.
//          Buttons and sliders should be prefixed with this string.
//
function Lab2(canvasID) {
  this.canvasID = canvasID;
  this.canvas = document.getElementById(canvasID);
  if (!this.canvas) {
    alert("Canvas ID '" + canvasID + "' not found.");
    return;
  }
  this.gl = WebGLUtils.setupWebGL(this.canvas);
  if (!this.gl) {
    alert("WebGL isn't available in this browser");
    return;
  }

  this.init();
}

// Define prototype values common to all Lab2 objects
Lab2.prototype.gl = null;

Lab2.prototype.toString = function () {
  return JSON.stringify(this);
};

Lab2.prototype.init = function () {
  var canvas = this.canvas;
  var gl = this.gl;
  var t = this;  // make available to event handlers

  // WebGL setup
  gl.viewport(0, 0, canvas.width, canvas.height);

  // Compile and link shaders
  this.shaderProgram = initShaders(gl, "vShader.glsl", "fShader.glsl");
  if (this.shaderProgram === null)
    return;
  gl.useProgram(this.shaderProgram);

  // Define names for colors
  var white = vec3(1.0, 1.0, 1.0);
  var red = vec3(1.0, 0.0, 0.0);
  var green = vec3(0.0, 1.0, 0.0);
  var blue = vec3(0.0, 0.0, 1.0);
  var yellow = vec3(1.0, 1.0, 0.0);
  
//  // Array of alternating initial vertex coordinates and colors for each vertex
//  this.vertexData = flatten([
//    vec3(0.0, 0.0, 0.0),
//    white,
//    vec3(0.5, 0.5, 0.0),
//    red,
//    vec3(-0.5, 0.5, 0.0),
//    green,
//    vec3(-0.5, -0.5, 0.0),
//    blue,
//    vec3(0.5, -0.5, 0.0),
//    yellow
//  ]);
  
  // Count of points in vertexData
  this.pointCount = 0;
  
  this.drawType;
  
  this.pArray = [];
  
  var floatSize = 4;  // size of gl.FLOAT in bytes
  // Load vertex data into WebGL buffer
  this.vertexCoordBuffer = gl.createBuffer();  // get unique buffer ID
  gl.bindBuffer(gl.ARRAY_BUFFER, this.vertexCoordBuffer);  // make this the active buffer
  this.gl.bufferData(gl.ARRAY_BUFFER, flatten(this.pArray), gl.STATIC_DRAW);  // write data to buffer

  // Define data layout in buffer for position.  Postions are 3 floats,
  // interleaved with 3 floats for colors, starting at beginning of buffer.
  this.vPosition = gl.getAttribLocation(this.shaderProgram, "vPosition");
  gl.vertexAttribPointer(this.vPosition, 3, gl.FLOAT, false, 6 * floatSize, 0);
  gl.enableVertexAttribArray(this.vPosition);

  // Define data layout in buffer for colors.  Colors are 3 floats,
  // interleaved with 3 floats for positions, starting after first position in buffer.
  this.vColor = gl.getAttribLocation(this.shaderProgram, "vColor");
  gl.vertexAttribPointer(this.vColor, 3, gl.FLOAT, false, 6 * floatSize, 3 * floatSize);
  gl.enableVertexAttribArray(this.vColor);

  // Define callback for change of slider value
  var sliderCallback = function (e) {
    // Update text display for slider
    var color = e.target.value;
    e.target.valueDisplay.textContent = color;

    // Re-render canvas
    requestAnimationFrame(render);
  };

  // Set up HTML user interface
  this.colors = ["r", "g", "b"];
  var rgbSliders = [];         // array of slider HTML elements
  var rgbSliderValues = [];    // array of slider value HTML elements

  // Set up an object with sliders for the three colors. The sliders are
  // accessed using "indices" of "r", "g", and "b".
  for (var i in this.colors) {
    var color = this.colors[i];
    var sliderID = this.canvasID + "-" + color + "-slider";
    rgbSliders[color] = document.getElementById(sliderID);
    if (rgbSliders[color] === null) {
      alert("Slider ID not found: " + sliderID);
      return;
    }
    var valueID = this.canvasID + "-" + color + "-value";
    rgbSliderValues[color] = document.getElementById(valueID);
    if (rgbSliders[color] === null) {
      alert("Slider value ID not found: " + sliderID);
      return;
    }
    rgbSliders[color].valueDisplay = rgbSliderValues[color];  // attach to slider
    
    // Set callback on slider input
    rgbSliders[color].addEventListener("input", sliderCallback);
  }
  this.rgbSliders = rgbSliders;

  var resetButton = document.getElementById(this.canvasID + "-reset-button");
  if (resetButton === null) {
    alert("Reset button ID not found: " + this.canvasID + "-reset-button");
    return;
  }

  // Set up callback to render a frame
  var render = function () {
    t.Render();
  };

  // Set up the callback for the reset button
  resetButton.addEventListener("click", function () {
    // Reset all the sliders to the middle value
    for (var i in rgbSliders) {
      rgbSliders[i].value = rgbSliders[i].max / 2.0;
      rgbSliders[i].valueDisplay.textContent =
              rgbSliders[i].valueAsNumber / rgbSliders[i].max;
    }
    t.pointCount = 0;
    t.pArray = [];
    requestAnimationFrame(render);
  });
  
  //Set up radio buttons
  var pointsButton = document.getElementById(this.canvasID + "-points");
  var linesButton = document.getElementById(this.canvasID + "-lines");
  var lsButton = document.getElementById(this.canvasID + "-linestrips");
  var trianglesButton = document.getElementById(this.canvasID + "-triangles");
  var tsButton = document.getElementById(this.canvasID + "-trianglestrips");
  var tfButton = document.getElementById(this.canvasID + "-trianglefan");

  // Set up mouse tracking
  var mouseX = document.getElementById(this.canvasID + "-mousex");
  var mouseY = document.getElementById(this.canvasID + "-mousey");
  //For Clip Coordinates
  var mouseClipX = document.getElementById(this.canvasID + "-mouseClipX");
  var mouseClipY = document.getElementById(this.canvasID + "-mouseClipY");
  //button
  var mouseButton = document.getElementById(this.canvasID + "-mousebutton");
  this.mouseDown = [ false, false, false ];  // track mouse button state
  mouseButton.textContent = this.mouseDown;
  if (mouseX === null || mouseY === null || mouseButton === null) {
    alert("Mouse output HTML IDs not found");
    return;
  }

  // Add mouse event handlers
  canvas.addEventListener("mousedown", function (e) {
    t.mouseDown[e.button] = true;
    mouseButton.textContent = t.mouseDown;
  });
  canvas.addEventListener("mouseup", function (e) {
    t.mouseDown[e.button] = false;
    mouseButton.textContent = t.mouseDown;
    
    //Push point location to array
    var x = ((2.0 * (e.pageX - e.target.offsetLeft)) / (canvas.width - 1)) - 1.0;;
    var y = 1.0 - ((2.0 * (e.pageY - e.target.offsetTop)) / (canvas.height -1));
    t.pArray.push(vec3(x, y, 0.0));
    t.pointCount++; //Increment point count
    var colorV = t.getSliderColor();
    t.pArray.push(colorV);
    gl.bufferData(gl.ARRAY_BUFFER, flatten(t.pArray), gl.STATIC_DRAW);  // write data to buffer
    // rendering
    requestAnimationFrame(render);
  });
  canvas.addEventListener("mousemove", function (e) {
    var x = e.pageX - e.target.offsetLeft;
    var y = e.pageY - e.target.offsetTop;
    mouseX.textContent = x;
    mouseY.textContent = y;
    //Clip coordinates
    var ClipX = ((2.0 * x) / (canvas.width - 1)) - 1.0;
    var ClipY = 1.0 - ((2.0 * y) / (canvas.height -1));
    mouseClipX.textContent = ClipX.toFixed(3);
    mouseClipY.textContent = ClipY.toFixed(3);
  });
};

/**
 * GetSliderColors - get the current RGB color represented by the sliders
 *   as a vec3.
 *   
 * @returns {vec3} current slider color
 */
Lab2.prototype.getSliderColor = function () {
  // Build an array of color values based on the current slider colors
  var colorValues = [];
  for (var i in this.colors) {
    var color = this.colors[i];
    var colorValue = this.rgbSliders[color].valueAsNumber;
    colorValues[i] = colorValue;
  }

  return vec3(colorValues);
};

/**
 * Render - draw the frame
 *
 */
Lab2.prototype.Render = function () {
  var gl = this.gl;
  gl.clearColor(0.0, 0.0, 0.0, 1.0);
  gl.clear(gl.COLOR_BUFFER_BIT);
if(document.getElementById("lines").checked)
    {
        gl.drawArrays(gl.LINES, 0, this.pointCount);
    }
    if(document.getElementById("linestrip").checked)
    {
        gl.drawArrays(gl.LINE_STRIP, 0, this.pointCount);
    }
    if(document.getElementById("lineloop").checked)
    {
        gl.drawArrays(gl.LINE_LOOP, 0, this.pointCount);
    }
    if(document.getElementById("triangle").checked)
    {
        gl.drawArrays(gl.TRIANGLES, 0, this.pointCount);
    }
    if(document.getElementById("trianglestrip").checked)
    {
        gl.drawArrays(gl.TRIANGLE_STRIP, 0, this.pointCount);
    }
    if(document.getElementById("trianglefan").checked)
    {
        gl.drawArrays(gl.TRIANGLE_FAN, 0, this.pointCount);
    }
     else
    {
        gl.drawArrays(gl.POINTS, 0, this.pointCount);
    }
};