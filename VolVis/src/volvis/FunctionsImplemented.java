/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package volvis;

import util.VectorMath;
import volume.VoxelGradient;

/**
 *
 * @author Administrator
 */
public class FunctionsImplemented
{
    // TODO 1:
    private short getVoxelTrilinear(double[] coord)
    {
        // TODO 1: Implement Tri-Linear interpolation and use it in your code
        // instead of getVoxel().
        
        // Each of the coordinate can take values only between
        // 0 and (max_value - 1)
        if (coord[0] < 0 || coord[0] > volume.getDimX() - 1 ||
            coord[1] < 0 || coord[1] > volume.getDimY() - 1 ||
            coord[2] < 0 || coord[2] > volume.getDimZ() - 1)
        {
            return 0;
        }
        
        // Further code is implemented by referring to the following article
        // https://en.wikipedia.org/wiki/Trilinear_interpolation
        
        //Calculate x0, y0, z0 and x1, y1, z1
        int x0 = (int) Math.floor(coord[0]);
        int y0 = (int) Math.floor(coord[1]);
        int z0 = (int) Math.floor(coord[2]);
        
        int x1 = (int) Math.ceil(coord[0]);
        int y1 = (int) Math.ceil(coord[1]);
        int z1 = (int) Math.ceil(coord[2]);
        
        //Calculating the co-efficients : alpha, beta and gamma
        double alpha = (coord[0] - x0) / (x1 - x0);
        double beta = (coord[1] - y0) / (y1 - y0);
        double gamma = (coord[2] - z0) / (z1 - z0);
        
        double c000 = volume.getVoxel(x0, y0, z0); // [0,0,0]
        double c001 = volume.getVoxel(x0, y0, z1); // [0,0,1]
        double c010 = volume.getVoxel(x0, y1, z0); // [0,1,0]
        double c011 = volume.getVoxel(x0, y1, z1); // [0,1,1]
        double c100 = volume.getVoxel(x1, y0, z0); // [1,0,0]
        double c101 = volume.getVoxel(x1, y0, z1); // [1,0,1]
        double c110 = volume.getVoxel(x1, y1, z0); // [1,1,0]
        double c111 = volume.getVoxel(x1, y1, z1); // [1,1,1]
        
        // Final computation of Tri-Linear Interpolation
        short interpolated_result = (short) Math.round(
                (1 - alpha) * (1 - beta) * (1 - gamma) * c000 +
                alpha * (1 - beta) * (1 - gamma) * c100 +
                (1 - alpha) * beta * (1 - gamma) * c010 +
                alpha * beta * (1 - gamma) * c110 +
                (1 - alpha) * (1 - beta) * gamma * c001 +
                alpha * (1 - beta) * gamma * c101 + 
                (1 - alpha) * beta * gamma * c011 + 
                alpha * beta * gamma * c111);

        return interpolated_result;
    }
    
    // TODO 2:
    private int traceRayComposite(double[] entryPoint, double[] exitPoint, double[] rayVector, double sampleStep, boolean frontBool)
    {
        // the light vector is directed toward the view point (which is the source of the light)
        // another light vector would be possible
        double[] lightVector = new double[3];
        VectorMath.setVector(lightVector, rayVector[0], rayVector[1], rayVector[2]);

        // creating and initializing the view vector
        double[] viewVector = new double[3];
        VectorMath.setVector(viewVector, -rayVector[0], -rayVector[1], -rayVector[2]);

        TFColor voxelColor = new TFColor(); // color before Phong shading
        TFColor voxelColor1 = new TFColor(); // color after Phong shading
        TFColor currentColor = new TFColor(); // color from TFunc2D
        TFColor colorAux = new TFColor(); // accumulated color
        
        // initialize colorAux
        colorAux.r = 0.0;
        colorAux.g = 0.0;
        colorAux.b = 0.0;
        colorAux.a = 0.0;

        //compute the increment and the number of samples (inverse: back-to-front)
        double[] increments = new double[3];
        VectorMath.setVector(increments, -rayVector[0] * sampleStep, -rayVector[1] * sampleStep, -rayVector[2] * sampleStep);

        // Compute the number of times we need to sample
        double distance = VectorMath.distance(entryPoint, exitPoint);
        int nrSamples = 1 + (int) Math.floor(distance / sampleStep);

        //the current position is initialized as the EXIT point
        double[] currentPos = new double[3];
        VectorMath.setVector(currentPos, exitPoint[0], exitPoint[1], exitPoint[2]);

        // TODO 2: To be Implemented this function. Now, it just gives back a constant color depending on the mode

        // Choosing Mode, TransferFunction, TransferFunction2D based on front/back side
        RaycastMode currMode = getMode(frontBool);
        TransferFunction currTFunc = getTFunc(frontBool);
        TransferFunction2D currTFunc2D = getTFunc2D(frontBool);

        do
        {
            // get the gradient and the intensity of the voxel
            int voxelValue = getVoxelTrilinear(currentPos);
            VoxelGradient gradient = getGradientTrilinear(currentPos);

            // 1D transfer function
            if (currMode == RaycastMode.COMPOSITING)
            {
                voxelColor = currTFunc.getColor(voxelValue);
            }

            if (currMode == RaycastMode.TRANSFER2D) 
            {
                currentColor = currTFunc2D.color;
                
                voxelColor.r = currentColor.r;
                voxelColor.g = currentColor.g;
                voxelColor.b = currentColor.b;
                
                // calling the compueOpacity2DTF function
                voxelColor.a = currentColor.a * computeOpacity2DTF(currTFunc2D.baseIntensity, 
                                                  currTFunc2D.radius,
                                                  voxelValue,
                                                  gradient.mag);
            }

            // If Phong Shading enabled
            if (shadingMode)
            {
                // Shading mode on
                //calculate the color and opacity after the application of phong shading
                currentColor = computePhongShading(voxelColor, gradient, lightVector, rayVector);

                //assign the color to our variable voxelColor1
                voxelColor1.r = currentColor.r;
                voxelColor1.g = currentColor.g;
                voxelColor1.b = currentColor.b;
                voxelColor1.a = currentColor.a;
            } else {
                voxelColor1.r = voxelColor.r;
                voxelColor1.g = voxelColor.g;
                voxelColor1.b = voxelColor.b;
                voxelColor1.a = voxelColor.a;
            }
            
            // compositing iteration, update accumulated color
            colorAux = compositingFromBackToFront(colorAux, voxelColor1);

            for (int i = 0; i < 3; i++)
            {
                currentPos[i] += increments[i];
            }

            nrSamples--;
        } while (nrSamples > 0);

        //computes the color
        int color = computePackedPixelColor(colorAux.r, colorAux.g, 
                                            colorAux.b, colorAux.a);
        return color;
    }
    
    private TFColor compositingFromBackToFront(TFColor compositedColor, TFColor currentColor)
    {
        double alpha = currentColor.a;
        compositedColor.r = currentColor.r * alpha + (1.0 - alpha) * compositedColor.r;
        compositedColor.g = currentColor.g * alpha + (1.0 - alpha) * compositedColor.g;
        compositedColor.b = currentColor.b * alpha + (1.0 - alpha) * compositedColor.b;
        compositedColor.a = alpha + (1.0 - alpha) * compositedColor.a;

        return compositedColor;
   }
   
    private TFColor compositingFromBackToFront1(TFColor compositedColor, TFColor currentColor, double alpha)
    {
        compositedColor.r = currentColor.r * alpha + (1.0 - alpha) * compositedColor.r;
        compositedColor.g = currentColor.g * alpha + (1.0 - alpha) * compositedColor.g;
        compositedColor.b = currentColor.b * alpha + (1.0 - alpha) * compositedColor.b;
        compositedColor.a = alpha + (1.0 - alpha) * compositedColor.a;
        
        return compositedColor;
   }
    
    
    // TODO 3:
    private int traceRayIso(double[] entryPoint, double[] exitPoint, double[] rayVector, double sampleStep, boolean frontBool)
    {

        double[] lightVector = new double[3];
        //We define the light vector as directed toward the view point (which is the source of the light)
        // another light vector would be possible
        
        VectorMath.setVector(lightVector, rayVector[0], rayVector[1], rayVector[2]);
        
        // TODO 3: Implement isosurface rendering.
        //Initialization of the colors as floating point values
        double r, g, b;
        r = g = b = 0.0;
        double alpha = 0.0;
        double opacity = 0;

        //Now we calculate the increase in samples
        double[] increment = new double[3];
        VectorMath.setVector(increment, - rayVector[0] * sampleStep, - rayVector[1] * sampleStep, - rayVector[2] * sampleStep);
                
        // Compute the number of times we need to sample
        double distance = VectorMath.distance(entryPoint, exitPoint);
        int nrSamples = 1 + (int) Math.floor(VectorMath.distance(entryPoint, exitPoint) / sampleStep);
                
        //the current position is the exit point
        double[] currentPos = new double[3];
        VectorMath.setVector(currentPos, exitPoint[0], exitPoint[1], exitPoint[2]);
        
        // Extract isoValue and color
        float isoValue = getIsoValue(frontBool);
        TFColor isoColor = getIsoColor(frontBool);
        
        do {
			// calling tri-linear interpolation method to get the interpolated voxel
            int value = getVoxelTrilinear(currentPos);

            if (value > isoValue) 
			{    
                r = isoColor.r;
                g = isoColor.g;
                b = isoColor.b;
                alpha = isoColor.a;
			
				// check if Shading should be applied
				if (shadingMode) 
				{
					TFColor d = new TFColor(r,g,b,alpha);
					VoxelGradient gradient = getGradientTrilinear(currentPos);
				
					TFColor new_color = computePhongShading(d, gradient, lightVector, rayVector);
					r = new_color.r;
					g = new_color.g;
					b = new_color.b;
					alpha = new_color.a;
				}
            }
          
            for (int j = 0; j < 3; j++) 
			{
                currentPos[j] += increment[j];
            }
           
           nrSamples--;
        } while (nrSamples > 0);
      
        
        //computes the color
        int color = computePackedPixelColor(r, g, b, alpha);
        return color;
    }
    
    // TODO 4:
    private void compute() 
    {

        for (int i=0; i<data.length; i++)
        {
            data[i] = zero;
        }
       
        for (int z=1; z<dimZ-1; z++)
        {
            for (int y=1; y<dimY-1; y++)
            {
                for (int x=1; x<dimX-1; x++)
                {
                    float x0 = (volume.getVoxel(x+1, y, z) - volume.getVoxel(x-1, y, z))/2.0f;
                    float y0 = (volume.getVoxel(x, y+1, z) - volume.getVoxel(x, y-1, z))/2.0f;
                    float z0 = (volume.getVoxel(x, y, z+1) - volume.getVoxel(x, y, z-1))/2.0f;
                    VoxelGradient newvalue = new VoxelGradient(x0, y0, z0);
                    setGradient(x, y, z, newvalue);
                }
            }
        }
     
     }
    
    // TODO 5:
    
    
    
    // TODO 6:
    private VoxelGradient getGradientTrilinear(double[] coord) 
    {
        // TODO 6: Implement Tri-linear interpolation for gradients
        
        // Get the coordinates
        double dx = coord[0], dy = coord[1], dz = coord[2];

        // Verify they are inside the volume gradient
        if (dx < 0 || dx > (gradients.getDimX() - 2) ||
            dy < 0 || dy > (gradients.getDimY() - 2) ||
            dz < 0 || dz > (gradients.getDimZ() - 2))
        {
            // If not, just return a zero gradient
            return ZERO_GRADIENT;
        }
        
        //Calculate x0, y0, z0 and x1, y1, z1
        int x0 = (int) Math.floor(dx);
        int y0 = (int) Math.floor(dy);
        int z0 = (int) Math.floor(dz);
        
        int x1 = (int) Math.ceil(dx);
        int y1 = (int) Math.ceil(dy);
        int z1 = (int) Math.ceil(dz);
        
        //Calculating the co-efficients : alpha, beta and gamma
        float alpha = (float)(dx - x0) / (x1 - x0);
        float beta = (float)(dy - y0) / (y1 - y0);
        float gamma = (float)(dz - z0) / (z1 - z0);
        
        // Computing the tri-linear interpolation of the X-component of the Gradient
        float gradientX = 
                (1 - alpha) * (1 - beta) * (1 - gamma) * gradients.getGradient(x0, y0, z0).x +
                alpha * (1 - beta) * (1 - gamma) * gradients.getGradient(x1, y0, z0).x +
                (1 - alpha) * beta * (1 - gamma) * gradients.getGradient(x0, y1, z0).x +
                alpha * beta * (1 - gamma) * gradients.getGradient(x1, y1, z0).x +
                (1 - alpha) * (1 - beta) * gamma * gradients.getGradient(x0, y0, z1).x +
                alpha * (1 - beta) * gamma * gradients.getGradient(x1, y0, z1).x + 
                (1 - alpha) * beta * gamma * gradients.getGradient(x0, y1, z1).x + 
                alpha * beta * gamma * gradients.getGradient(x1, y1, z1).x;
        
        // Computing the tri-linear interpolation of the Y-component of the Gradient
        float gradientY = 
                (1 - alpha) * (1 - beta) * (1 - gamma) * gradients.getGradient(x0, y0, z0).y +
                alpha * (1 - beta) * (1 - gamma) * gradients.getGradient(x1, y0, z0).y +
                (1 - alpha) * beta * (1 - gamma) * gradients.getGradient(x0, y1, z0).y +
                alpha * beta * (1 - gamma) * gradients.getGradient(x1, y1, z0).y +
                (1 - alpha) * (1 - beta) * gamma * gradients.getGradient(x0, y0, z1).y +
                alpha * (1 - beta) * gamma * gradients.getGradient(x1, y0, z1).y + 
                (1 - alpha) * beta * gamma * gradients.getGradient(x0, y1, z1).y + 
                alpha * beta * gamma * gradients.getGradient(x1, y1, z1).y;
        
        // Computing the tri-linear interpolation of the Z-component of the Gradient
        float gradientZ = 
                (1 - alpha) * (1 - beta) * (1 - gamma) * gradients.getGradient(x0, y0, z0).z +
                alpha * (1 - beta) * (1 - gamma) * gradients.getGradient(x1, y0, z0).z +
                (1 - alpha) * beta * (1 - gamma) * gradients.getGradient(x0, y1, z0).z +
                alpha * beta * (1 - gamma) * gradients.getGradient(x1, y1, z0).z +
                (1 - alpha) * (1 - beta) * gamma * gradients.getGradient(x0, y0, z1).z +
                alpha * (1 - beta) * gamma * gradients.getGradient(x1, y0, z1).z + 
                (1 - alpha) * beta * gamma * gradients.getGradient(x0, y1, z1).z + 
                alpha * beta * gamma * gradients.getGradient(x1, y1, z1).z;
        
        VoxelGradient resultantGradient = new VoxelGradient(gradientX,
                                                            gradientY,
                                                            gradientZ);
        
        return resultantGradient;
    }
    
    // TODO 7:
    private TFColor computePhongShading(TFColor voxel_color, VoxelGradient gradient,
            double[] lightVector, double[] rayVector)
    {

        // TODO 7: Implement Phong Shading.

        // Return the exisitng voxel color if the voxel color is zero
        if(gradient.mag == 0)
        {
            return voxel_color;
        }

        double ka = 0.1; // ambient
        double kd = 0.7; // diffuse
        double ks = 0.2; // specular

        int alpha = 100;

        double[] lightProperty_a = {1,1,1};
        double[] lightProperty_d = {1,1,1};
        double[] lightProperty_s = {1,1,1};

        TFColor color;

        //fish doesn't exist
        //if (gradient.x == 0 && gradient.y == 0 && gradient.z == 0) return voxel_color;
        if(gradient.mag == 0)   return voxel_color;

        //norm gradiente
        double[] gradVec = {gradient.x / gradient.mag, gradient.y / gradient.mag, gradient.z / gradient.mag};

        //norm light
        double lightNorm = VectorMath.length(lightVector);  //Math.sqrt(Math.pow(lightVector[0], 2) + Math.pow(lightVector[1], 2) + Math.pow(lightVector[2], 2));
        if(lightNorm == 0) return voxel_color;
        double[] lightNormVector = {lightVector[0] / lightNorm, lightVector[1] / lightNorm, lightVector[2] / lightNorm};

        //norm ray
        double rayNorm = VectorMath.length(rayVector);  //Math.sqrt(Math.pow(rayVector[0], 2) + Math.pow(rayVector[1], 2) + Math.pow(rayVector[2], 2));
        if(rayNorm == 0) return voxel_color;
        double[] rayNormVector = {rayVector[0] / rayNorm, rayVector[1] / rayNorm, rayVector[2] / rayNorm};

        //cos 1
        double cos1 = VectorMath.dotproduct(lightNormVector, gradVec) > 0 ? VectorMath.dotproduct(lightNormVector, gradVec) : 0;

        //R
        double[] twice_gradVec = {gradVec[0] * 2, gradVec[1] * 2, gradVec[2] * 2};
        double x = VectorMath.dotproduct(twice_gradVec, lightNormVector);
        double[] x_gradVec = {gradVec[0] * x, gradVec[1] * x, gradVec[2] * x};
        double[] R = {x_gradVec[0] - lightNormVector[0], x_gradVec[1] - lightNormVector[1], x_gradVec[2] - lightNormVector[2]};
        double cos2 = VectorMath.dotproduct(rayNormVector, R) > 0 ? VectorMath.dotproduct(rayNormVector, R): 0;

        double r = lightProperty_a[0] * ka * voxel_color.r +
                   lightProperty_d[0] * kd * voxel_color.r * cos1 +
                   lightProperty_s[0] * ks * voxel_color.r * Math.pow(cos2, alpha);
        double g = lightProperty_a[1] * ka * voxel_color.g +
                   lightProperty_d[1] * kd * voxel_color.g * cos1 +
                   lightProperty_s[1] * ks * voxel_color.g * Math.pow(cos2, alpha);
        double b = lightProperty_a[2] * ka * voxel_color.b +
                   lightProperty_d[2] * kd * voxel_color.b * cos1 +
                   lightProperty_s[2] * ks * voxel_color.b * Math.pow(cos2, alpha);

        if (r < 0) {
            r = 0;
        }
        if (g < 0) {
            g = 0;
        }
        if (b < 0) {
            b = 0;
        }
        if (r > 1) {
            r = 1;
        }
        if (g > 1) {
            g = 1;
        }
        if (b > 1) {
            b = 1;
        }

        color = new TFColor(r, g, b, voxel_color.a);
        return color;

    }
    
    // TODO 8:
    public double computeOpacity2DTF(double material_value, double material_r,
            double voxelValue, double gradMagnitude)
    {
        double opacity = 0.0;

        // adding the suggestion announced by Humberto
        double radius = material_r/gradients.getMaxGradientMagnitude();

        // TODO 8: Implement weight based opacity.

        // define the upper and lower threshold bounds for material_value

        double lowerThresholdForMaterialValue = voxelValue - radius * gradMagnitude;
        double upperThresholdForMaterialValue = voxelValue + radius * gradMagnitude;

        // compute the value of opacity differently for 3 different inequalities as mentioned in Levoy's paper

        if(material_value == voxelValue && gradMagnitude == 0)
        {
            opacity = 1.0;
        }
        else if(gradMagnitude > 0 &&
		lowerThresholdForMaterialValue <= material_value &&
		material_value <= upperThresholdForMaterialValue)
        {
            opacity = 1 - (1/radius) * Math.abs((voxelValue - material_value)/gradMagnitude);
        }
        else
        {
            opacity = 0.0;
        }

        return opacity;
    }
    
    // TODO 9:
    
    
}
