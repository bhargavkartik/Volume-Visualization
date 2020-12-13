package volvis;

import com.jogamp.opengl.GL;
import com.jogamp.opengl.GL2;
import com.jogamp.opengl.util.texture.Texture;
import com.jogamp.opengl.util.texture.awt.AWTTextureIO;
import gui.RaycastRendererPanel;
import gui.TransferFunction2DEditor;
import gui.TransferFunctionEditor;
import java.awt.image.BufferedImage;

import util.TFChangeListener;
import util.VectorMath;
import volume.GradientVolume;
import volume.Volume;
import volume.VoxelGradient;

/**
 * Raycast Renderer.
 *
 * @author Michel Westenberg
 * @author Anna Vilanova
 * @author Nicola Pezzotti
 * @author Humberto Garcia
 */
public class RaycastRenderer extends Renderer implements TFChangeListener {

    /**
     * Volume that is loaded and visualized.
     */
    private Volume volume = null;

    /**
     * Rendered image.
     */
    private BufferedImage image;

    /**
     * Gradient information of the loaded volume.
     */
    private GradientVolume gradients = null;

    /**
     * Reference to the GUI panel.
     */
    RaycastRendererPanel panelFront;

    /**
     * Transfer Function.
     */
    TransferFunction tFuncFront;

    /**
     * Reference to the GUI transfer function editor.
     */
    TransferFunctionEditor tfEditor;

    /**
     * Transfer Function 2D.
     */
    TransferFunction2D tFunc2DFront;

    /**
     * Reference to the GUI 2D transfer function editor.
     */
    TransferFunction2DEditor tfEditor2DFront;

    /**
     * Mode of our raycast. See {@link RaycastMode}
     */
    private RaycastMode modeFront;

    /**
     * Whether we are in cutting plane mode or not.
     */
    private boolean cuttingPlaneMode = false;

    /**
     * Whether we are in shading mode or not.
     */
    private boolean shadingMode = false;

    /**
     * Iso value to use in Isosurface rendering.
     */
    private float isoValueFront = 95f;

    /**
     * Color used for the isosurface rendering.
     */
    private TFColor isoColorFront;

    // Below cutting plane specific attributes
    /**
     * Cutting plane normal vector.
     */
    private final double[] planeNorm = new double[]{0d, 0d, 1d};

    /**
     * Cutting plane point.
     */
    private final double[] planePoint = new double[]{0d, 0d, 0d};

    /**
     * Back mode of our raycast for cutting plane.
     */
    private RaycastMode modeBack;

    /**
     * Iso value to use in Isosurface rendering for cutting plane.
     */
    private float isoValueBack = 95f;

    /**
     * Color used for the isosurface rendering for cutting plane.
     */
    private TFColor isoColorBack;

    /**
     * Transfer Function for cutting plane.
     */
    TransferFunction tFuncBack;

    /**
     * Reference to the GUI transfer function editor for cutting plane.
     */
    TransferFunctionEditor tfEditorBack;

    /**
     * Transfer Function 2D for cutting plane.
     */
    TransferFunction2D tFunc2DBack;

    /**
     * Reference to the GUI 2D transfer function editor for cutting plane.
     */
    TransferFunction2DEditor tfEditor2DBack;

    /**
     * Constant Zero gradient.
     */
    private final static VoxelGradient ZERO_GRADIENT = new VoxelGradient();

    /**
     * Gets the corresponding voxel using Nearest Neighbors.
     *
     * @param coord Pixel coordinate in 3D space of the voxel we want to get.
     * @return The voxel value.
     */
    private short getVoxel(double[] coord) {
        // Get coordinates
        double dx = coord[0], dy = coord[1], dz = coord[2];

        // Verify they are inside the volume
        if (dx < 0 || dx >= volume.getDimX() || dy < 0 || dy >= volume.getDimY()
                || dz < 0 || dz >= volume.getDimZ()) {

            // If not, jus return 0
            return 0;
        }

        // Get the closest x, y, z to dx, dy, dz that are integers
        // This is important as our data is discrete (not continuous)
        int x = (int) Math.floor(dx);
        int y = (int) Math.floor(dy);
        int z = (int) Math.floor(dz);

        // Finally, get the voxel from the Volume for the corresponding coordinates
        return volume.getVoxel(x, y, z);
    }

    /**
     * Gets the corresponding voxel using Tri-linear Interpolation.
     *
     * @param coord Pixel coordinate in 3D space of the voxel we want to get.
     * @return The voxel value.
     */
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

    /**
     * Gets the corresponding VoxelGradient using Nearest Neighbors.
     *
     * @param coord Pixel coordinate in 3D space of the voxel we want to get.
     * @return The voxel gradient.
     */
    private VoxelGradient getGradient(double[] coord) {
        // Get the coordinates
        double dx = coord[0], dy = coord[1], dz = coord[2];

        // Verify they are inside the volume gradient
        if (dx < 0 || dx > (gradients.getDimX() - 2) || dy < 0 || dy > (gradients.getDimY() - 2)
                || dz < 0 || dz > (gradients.getDimZ() - 2)) {

            // If not, just return a zero gradient
            return ZERO_GRADIENT;
        }

        // Get the closest x, y, z to dx, dy, dz that are integers
        // This is important as our data is discrete (not continuous)
        int x = (int) Math.round(dx);
        int y = (int) Math.round(dy);
        int z = (int) Math.round(dz);

        // Finally, get the gradient from GradientVolume for the corresponding coordinates
        return gradients.getGradient(x, y, z);
    }

    /**
     * Gets the corresponding VoxelGradient using Tri-linear interpolation.
     *
     * @param coord Pixel coordinate in 3D space of the voxel we want to get.
     * @return The voxel gradient.
     */
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

    /**
     * Updates {@link #image} attribute (result of rendering) using the slicing
     * technique.
     *
     * @param viewMatrix OpenGL View matrix {
     * @see
     * <a href="www.songho.ca/opengl/gl_transform.html#modelview">link</a>}.
     */
    private void slicer(double[] viewMatrix) {

        // Clear the image
        resetImage();

        // vector uVec and vVec define a plane through the origin,
        // perpendicular to the view vector viewVec which is going from the view point towards the object
        // uVec contains the up vector of the camera in world coordinates (image vertical)
        // vVec contains the horizontal vector in world coordinates (image horizontal)
        double[] viewVec = new double[3];
        double[] uVec = new double[3];
        double[] vVec = new double[3];
        VectorMath.setVector(viewVec, viewMatrix[2], viewMatrix[6], viewMatrix[10]);
        VectorMath.setVector(uVec, viewMatrix[0], viewMatrix[4], viewMatrix[8]);
        VectorMath.setVector(vVec, viewMatrix[1], viewMatrix[5], viewMatrix[9]);

        // We get the size of the image/texture we will be puting the result of the 
        // volume rendering operation.
        int imageW = image.getWidth();
        int imageH = image.getHeight();

        int[] imageCenter = new int[2];
        // Center of the image/texture 
        imageCenter[0] = imageW / 2;
        imageCenter[1] = imageH / 2;

        double[] pixelCoord = new double[3];
        double[] volumeCenter = new double[3];
        VectorMath.setVector(volumeCenter, volume.getDimX() / 2, volume.getDimY() / 2, volume.getDimZ() / 2);

        // sample on a plane through the origin of the volume data
        double max = volume.getMaximum();

        TFColor pixelColor = new TFColor();
        // Auxiliar color
        TFColor colorAux;

        for (int j = imageCenter[1] - imageH / 2; j < imageCenter[1] + imageH / 2; j++) {
            for (int i = imageCenter[0] - imageW / 2; i < imageCenter[0] + imageW / 2; i++) {
                // computes the pixelCoord which contains the 3D coordinates of the pixels (i,j)
                computePixelCoordinatesFloat(pixelCoord, volumeCenter, uVec, vVec, i, j);

                //int val = getVoxel(pixelCoord);
                //NOTE: you have to implement this function to get the tri-linear interpolation
                int val = getVoxelTrilinear(pixelCoord);

                // Map the intensity to a grey value by linear scaling
                pixelColor.r = val / max;
                pixelColor.g = pixelColor.r;
                pixelColor.b = pixelColor.r;
                pixelColor.a = val > 0 ? 1.0 : 0.0;  // this makes intensity 0 completely transparent and the rest opaque
                // Alternatively, apply the transfer function to obtain a color
                // pixelColor = tFuncFront.getColor(val);

                //BufferedImage/image/texture expects a pixel color packed as ARGB in an int
                //use the function computeImageColor to convert your double color in the range 0-1 to the format need by the image
                int packedPixelColor = computePackedPixelColor(pixelColor.r, pixelColor.g, pixelColor.b, pixelColor.a);
                image.setRGB(i, j, packedPixelColor);
            }
        }
    }

    /**
     * Do NOT modify this function.
     *
     * Updates {@link #image} attribute (result of rendering) using MIP
     * raycasting. It returns the color assigned to a ray/pixel given its
     * starting and ending points, and the direction of the ray.
     *
     * @param entryPoint Starting point of the ray.
     * @param exitPoint Last point of the ray.
     * @param rayVector Direction of the ray.
     * @param sampleStep Sample step of the ray.
     * @return Color assigned to a ray/pixel.
     */
    private int traceRayMIP(double[] entryPoint, double[] exitPoint, double[] rayVector, double sampleStep) {
        //compute the increment and the number of samples
        double[] increments = new double[3];
        VectorMath.setVector(increments, rayVector[0] * sampleStep, rayVector[1] * sampleStep, rayVector[2] * sampleStep);

        // Compute the number of times we need to sample
        double distance = VectorMath.distance(entryPoint, exitPoint);
        int nrSamples = 1 + (int) Math.floor(VectorMath.distance(entryPoint, exitPoint) / sampleStep);

        //the current position is initialized as the entry point
        double[] currentPos = new double[3];
        VectorMath.setVector(currentPos, entryPoint[0], entryPoint[1], entryPoint[2]);

        double maximum = 0;
        do {
            double value = getVoxelTrilinear(currentPos) / 255.;
            if (value > maximum) {
                maximum = value;
            }
            for (int i = 0; i < 3; i++) {
                currentPos[i] += increments[i];
            }
            nrSamples--;
        } while (nrSamples > 0);

        double alpha;
        double r, g, b;
        if (maximum > 0.0) { // if the maximum = 0 make the voxel transparent
            alpha = 1.0;
        } else {
            alpha = 0.0;
        }
        r = g = b = maximum;
        int color = computePackedPixelColor(r, g, b, alpha);
        return color;
    }
    
    // Utility functions for cutting plane
    public TFColor getIsoColor(boolean frontBool) {
        if (frontBool) {
            return isoColorFront;
        } else {
            return isoColorBack;
        }
    }

    public float getIsoValue(boolean frontBool) {
        if (frontBool) {
            return isoValueFront;
        } else {
            return isoValueBack;
        }
    }

    public RaycastMode getMode(boolean frontBool) {
        if (frontBool) {
            return modeFront;
        } else {
            return modeBack;
        }
    }

    public TransferFunction getTFunc(boolean frontBool) {
        if (frontBool) {
            return tFuncFront;
        } else {
            return tFuncBack;
        }
    }

    public TransferFunction2D getTFunc2D(boolean frontBool) {
        if (frontBool) {
            return tFunc2DFront;
        } else {
            return tFunc2DBack;
        }
    }


    /**
     *
     * Updates {@link #image} attribute (result of rendering) using the
     * Isosurface raycasting. It returns the color assigned to a ray/pixel given
     * its starting and ending points, and the direction of the ray.
     *
     * @param entryPoint Starting point of the ray.
     * @param exitPoint Last point of the ray.
     * @param rayVector Direction of the ray.
     * @param sampleStep Sample step of the ray.
     * @return Color assigned to a ray/pixel.
     */
    private int traceRayIso(double[] entryPoint, double[] exitPoint, double[] rayVector, double sampleStep, boolean frontBool) {

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
        VectorMath.setVector(increment, -rayVector[0] * sampleStep, -rayVector[1] * sampleStep, -rayVector[2] * sampleStep);

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

            if (value > isoValue) {
                r = isoColor.r;
                g = isoColor.g;
                b = isoColor.b;
                alpha = isoColor.a;

                // check if Shading should be applied
                if (shadingMode) {
                    TFColor d = new TFColor(r, g, b, alpha);
                    VoxelGradient gradient = getGradientTrilinear(currentPos);

                    TFColor new_color = computePhongShading(d, gradient, lightVector, rayVector);
                    r = new_color.r;
                    g = new_color.g;
                    b = new_color.b;
                    alpha = new_color.a;
                }
            }

            for (int j = 0; j < 3; j++) {
                currentPos[j] += increment[j];
            }

            nrSamples--;
        } while (nrSamples > 0);

        //computes the color
        int color = computePackedPixelColor(r, g, b, alpha);
        return color;
    }

    /**
     *
     * Updates {@link #image} attribute (result of rendering) using the
     * compositing/accumulated raycasting. It returns the color assigned to a
     * ray/pixel given its starting and ending points, and the direction of the
     * ray.
     *
     * Ray must be sampled with a distance defined by sampleStep.
     *
     * @param entryPoint Starting point of the ray.
     * @param exitPoint Last point of the ray.
     * @param rayVector Direction of the ray.
     * @param sampleStep Sample step of the ray.
     * @return Color assigned to a ray/pixel.
     */
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

    /*
     * Compute the back to front compositing
     * @param compositedColor The Composited Color
     * @param currentColor The Current Color
     * @param alpha Transparency
     * @param compositedColor Returns the Composited color
     */

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

    /**
     * Compute Phong Shading given the voxel color (material color), gradient,
     * light vector and view vector.
     *
     * @param voxel_color Voxel color (material color).
     * @param gradient Gradient voxel.
     * @param lightVec Light vector.
     * @param rayVector View vector.
     * @return Computed color for Phong Shading.
     */
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

        // if (gradient.x == 0 && gradient.y == 0 && gradient.z == 0) return voxel_color;
        if(gradient.mag == 0)   return voxel_color;

        // norm gradiente
        double[] gradVec = {gradient.x / gradient.mag, gradient.y / gradient.mag, gradient.z / gradient.mag};

        // norm light
        double lightNorm = VectorMath.length(lightVector);  //Math.sqrt(Math.pow(lightVector[0], 2) + Math.pow(lightVector[1], 2) + Math.pow(lightVector[2], 2));
        if(lightNorm == 0) return voxel_color;
        double[] lightNormVector = {lightVector[0] / lightNorm, lightVector[1] / lightNorm, lightVector[2] / lightNorm};

        // norm ray
        double rayNorm = VectorMath.length(rayVector);  //Math.sqrt(Math.pow(rayVector[0], 2) + Math.pow(rayVector[1], 2) + Math.pow(rayVector[2], 2));
        if(rayNorm == 0) return voxel_color;
        double[] rayNormVector = {rayVector[0] / rayNorm, rayVector[1] / rayNorm, rayVector[2] / rayNorm};

        // cos 1
        double cos1 = VectorMath.dotproduct(lightNormVector, gradVec) > 0 ? VectorMath.dotproduct(lightNormVector, gradVec) : 0;

        // R
        double[] twice_gradVec = {gradVec[0] * 2, gradVec[1] * 2, gradVec[2] * 2};
        double x = VectorMath.dotproduct(twice_gradVec, lightNormVector);
        double[] x_gradVec = {gradVec[0] * x, gradVec[1] * x, gradVec[2] * x};
        double[] R = {x_gradVec[0] - lightNormVector[0], x_gradVec[1] - lightNormVector[1], x_gradVec[2] - lightNormVector[2]};
        double cos2 = VectorMath.dotproduct(rayNormVector, R) > 0 ? VectorMath.dotproduct(rayNormVector, R): 0;

        double r = ka * voxel_color.r +
                   kd * voxel_color.r * cos1 +
                   lightProperty_s[0] * ks * Math.pow(cos2, alpha);
        double g = ka * voxel_color.g +
                   kd * voxel_color.g * cos1 +
                   lightProperty_s[1] * ks * Math.pow(cos2, alpha);
        double b = ka * voxel_color.b +
                   kd * voxel_color.b * cos1 +
                   lightProperty_s[2] * ks * Math.pow(cos2, alpha);
        
        // r, g, b should be between 0.0 to 1.0
        if (r < 0) { r = 0; }
        if (g < 0) { g = 0; }
        if (b < 0) { b = 0; }
        if (r > 1) { r = 1; }
        if (g > 1) { g = 1; }
        if (b > 1) { b = 1; }

        color = new TFColor(r, g, b, voxel_color.a);
        return color;

    }

    /**
     * Implements the basic tracing of rays through the image given the camera
     * transformation. It calls the functions depending on the raycasting mode.
     *
     * @param viewMatrix
     */
    void raycast(double[] viewMatrix)
    {
        //data allocation
        double[] viewVec = new double[3];
        double[] uVec = new double[3];
        double[] vVec = new double[3];
        double[] pixelCoord = new double[3];
        double[] entryPoint = new double[3];
        double[] exitPoint = new double[3];

        // TODO 5: Limited modification is needed
        // increment in the pixel domain in pixel units
        int increment;
        // sample step in voxel units
        int sampleStep;

        //System.out.println(interactiveMode);
        if (interactiveMode) {
            increment = 2;
            sampleStep = 3;
        } else {
            increment = 1;
            sampleStep = 1;
        }

        // reset the image to black
        resetImage();

        // vector uVec and vVec define a plane through the origin,
        // perpendicular to the view vector viewVec which is going from the view point towards the object
        // uVec contains the up vector of the camera in world coordinates (image vertical)
        // vVec contains the horizontal vector in world coordinates (image horizontal)
        VectorMath.setVector(viewVec, viewMatrix[2], viewMatrix[6], viewMatrix[10]);
        VectorMath.setVector(uVec, viewMatrix[0], viewMatrix[4], viewMatrix[8]);
        VectorMath.setVector(vVec, viewMatrix[1], viewMatrix[5], viewMatrix[9]);

        // We get the size of the image/texture we will be puting the result of the 
        // volume rendering operation.
        int imageW = image.getWidth();
        int imageH = image.getHeight();

        int[] imageCenter = new int[2];
        // Center of the image/texture 
        imageCenter[0] = imageW / 2;
        imageCenter[1] = imageH / 2;

        //The rayVector is pointing towards the scene
        double[] rayVector = new double[3];
        rayVector[0] = -viewVec[0];
        rayVector[1] = -viewVec[1];
        rayVector[2] = -viewVec[2];

        // compute the volume center
        double[] volumeCenter = new double[3];
        VectorMath.setVector(volumeCenter, volume.getDimX() / 2, volume.getDimY() / 2, volume.getDimZ() / 2);
        
        // Vector from entru poin to the planePoint
        double[] diffVec = new double[3];
        
        // Whether the ray is for upper half-plane
        boolean frontBool = true;

        // ray computation for each pixel
        for (int j = imageCenter[1] - imageH / 2; j < imageCenter[1] + imageH / 2; j += increment) {
            for (int i = imageCenter[0] - imageW / 2; i < imageCenter[0] + imageW / 2; i += increment) {
                // compute starting points of rays in a plane shifted backwards to a position behind the data set
                computePixelCoordinatesBehindFloat(pixelCoord, viewVec, uVec, vVec, i, j);
                // compute the entry and exit point of the ray
                computeEntryAndExit(pixelCoord, rayVector, entryPoint, exitPoint);

                // TODO 9: Implement logic for cutting plane.
                if (this.isCuttingPlaneMode()) {
                    double dot = util.VectorMath.dotproduct(util.VectorMath.difference(entryPoint, planePoint, diffVec), planeNorm);
 
                    if (dot >= 0) {
                        frontBool = true;
                    } else {
                        frontBool = false;
                    }
                }
                RaycastMode currMode = getMode(frontBool);

                if ((entryPoint[0] > -1.0) && (exitPoint[0] > -1.0)) {
                    int val = 0;
                    switch (currMode) {
                        case COMPOSITING:
                        case TRANSFER2D:
                            val = traceRayComposite(entryPoint, exitPoint, rayVector, sampleStep, frontBool);
                            break;
                        case MIP:
                            val = traceRayMIP(entryPoint, exitPoint, rayVector, sampleStep);
                            break;
                        case ISO_SURFACE:
                            val = traceRayIso(entryPoint, exitPoint, rayVector, sampleStep, frontBool);
                            break;
                    }
                    for (int ii = i; ii < i + increment; ii++) {
                        for (int jj = j; jj < j + increment; jj++) {
                            image.setRGB(ii, jj, val);
                        }
                    }
                }

            }
        }
    }

    /**
     * Computes the opacity based on the value of the pixel and values of the
     * triangle widget. {@link #tFunc2DFront} contains the values of the base
     * intensity and radius. {@link TransferFunction2D#baseIntensity} and
     * {@link TransferFunction2D#radius} are in image intensity units.
     *
     * @param material_value Value of the material.
     * @param material_r Radius of the material.
     * @param voxelValue Voxel value.
     * @param gradMagnitude Gradient magnitude.
     * @return
     */
    public double computeOpacity2DTF(double material_value, double material_r,
            double voxelValue, double gradMagnitude)
    {
        double opacity = 0.0;

        // adding the suggestion announced by Humberto
        double radius = material_r/gradients.getMaxGradientMagnitude();

        // TODO 8: Implement weight based opacity.

        // TFColor baseColor=tFunc2DFront.color;

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

    /**
     * Class constructor. Initializes attributes.
     */
    public RaycastRenderer() {
        panelFront = new RaycastRendererPanel(this);
        panelFront.setSpeedLabel("0");

        isoColorFront = new TFColor();
        isoColorFront.r = 1.0;
        isoColorFront.g = 1.0;
        isoColorFront.b = 0.0;
        isoColorFront.a = 1.0;

        isoColorBack = new TFColor();
        isoColorBack.r = 1.0;
        isoColorBack.g = 1.0;
        isoColorBack.b = 0.0;
        isoColorBack.a = 1.0;

        modeFront = RaycastMode.SLICER;
        modeBack = RaycastMode.SLICER;
    }

    /**
     * Sets the volume to be visualized. It creates the Image buffer for the
     * size of the volume. Initializes the transfers functions
     *
     * @param vol Volume to be visualized.
     */
    public void setVolume(Volume vol) {
        System.out.println("Assigning volume");
        volume = vol;

        System.out.println("Computing gradients");
        gradients = new GradientVolume(vol);

        // set up image for storing the resulting rendering
        // the image width and height are equal to the length of the volume diagonal
        int imageSize = (int) Math.floor(Math.sqrt(vol.getDimX() * vol.getDimX() + vol.getDimY() * vol.getDimY()
                + vol.getDimZ() * vol.getDimZ()));
        if (imageSize % 2 != 0) {
            imageSize = imageSize + 1;
        }

        image = new BufferedImage(imageSize, imageSize, BufferedImage.TYPE_INT_ARGB);

        // Initialize transfer function and GUI panels
        tFuncFront = new TransferFunction(volume.getMinimum(), volume.getMaximum());
        tFuncFront.setTestFunc();
        tFuncFront.addTFChangeListener(this);
        tfEditor = new TransferFunctionEditor(tFuncFront, volume.getHistogram());

        tFunc2DFront = new TransferFunction2D((short) (volume.getMaximum() / 2), 0.2 * volume.getMaximum());
        tfEditor2DFront = new TransferFunction2DEditor(tFunc2DFront, volume, gradients);
        tfEditor2DFront.addTFChangeListener(this);

        // Initialize transfer function and GUI panels for cutting plane
        tFuncBack = new TransferFunction(volume.getMinimum(), volume.getMaximum());
        tFuncBack.setTestFunc();
        tFuncBack.addTFChangeListener(this);
        tfEditorBack = new TransferFunctionEditor(tFuncBack, volume.getHistogram());

        tFunc2DBack = new TransferFunction2D((short) (volume.getMaximum() / 2), 0.2 * volume.getMaximum());
        tfEditor2DBack = new TransferFunction2DEditor(tFunc2DBack, volume, gradients);
        tfEditor2DBack.addTFChangeListener(this);

        // Set plane point
        VectorMath.setVector(planePoint, volume.getDimX() / 2, volume.getDimY() / 2, volume.getDimZ() / 2);

        System.out.println("Finished initialization of RaycastRenderer");
    }

    /**
     * Do NOT modify.
     *
     * Visualizes the volume. It calls the corresponding render functions.
     *
     * @param gl OpenGL API.
     */
    @Override
    public void visualize(GL2 gl) {
        if (volume == null) {
            return;
        }

        drawBoundingBox(gl);

        // If mode is Cutting Plane, draw the cutting plane.
        if (cuttingPlaneMode) {
            drawCuttingPlane(gl);
        }

        gl.glGetDoublev(GL2.GL_MODELVIEW_MATRIX, _viewMatrix, 0);

        long startTime = System.currentTimeMillis();

        switch (modeFront) {
            case SLICER:
                slicer(_viewMatrix);
                break;
            default:
                // Default case raycast
                raycast(_viewMatrix);
                break;
        }

        long endTime = System.currentTimeMillis();
        double runningTime = (endTime - startTime);
        panelFront.setSpeedLabel(Double.toString(runningTime));

        Texture texture = AWTTextureIO.newTexture(gl.getGLProfile(), image, false);

        gl.glPushAttrib(GL2.GL_LIGHTING_BIT);
        gl.glDisable(GL2.GL_LIGHTING);
        gl.glEnable(GL.GL_BLEND);
        gl.glBlendFunc(GL.GL_SRC_ALPHA, GL.GL_ONE_MINUS_SRC_ALPHA);

        // draw rendered image as a billboard texture
        texture.enable(gl);
        texture.bind(gl);
        double halfWidth = image.getWidth() / 2.0;
        gl.glPushMatrix();
        gl.glLoadIdentity();
        gl.glBegin(GL2.GL_QUADS);
        gl.glColor4f(1.0f, 1.0f, 1.0f, 1.0f);
        gl.glTexCoord2d(0.0, 0.0);
        gl.glVertex3d(-halfWidth, -halfWidth, 0.0);
        gl.glTexCoord2d(0.0, 1.0);
        gl.glVertex3d(-halfWidth, halfWidth, 0.0);
        gl.glTexCoord2d(1.0, 1.0);
        gl.glVertex3d(halfWidth, halfWidth, 0.0);
        gl.glTexCoord2d(1.0, 0.0);
        gl.glVertex3d(halfWidth, -halfWidth, 0.0);
        gl.glEnd();
        texture.disable(gl);
        texture.destroy(gl);
        gl.glPopMatrix();

        gl.glPopAttrib();

        if (gl.glGetError() > 0) {
            System.out.println("some OpenGL error: " + gl.glGetError());
        }
    }

    public RaycastMode getRaycastMode() {
        return modeFront;
    }

    /**
     * Sets the raycast mode to the specified one.
     *
     * @param mode New Raycast mode.
     */
    public void setRaycastModeFront(RaycastMode mode) {
        this.modeFront = mode;
    }

    public void setRaycastModeBack(RaycastMode mode) {
        this.modeBack = mode;
    }

    @Override
    public void changed() {
        for (int i = 0; i < listeners.size(); i++) {
            listeners.get(i).changed();
        }
    }

    /**
     * Do NOT modify.
     *
     * Updates the vectors that represent the cutting plane.
     *
     * @param d View Matrix.
     */
    public void updateCuttingPlaneVectors(double[] d) {
        VectorMath.setVector(_planeU, d[1], d[5], d[9]);
        VectorMath.setVector(_planeV, d[2], d[6], d[10]);
        VectorMath.setVector(planeNorm, d[0], d[4], d[8]);
    }

    /**
     * Sets the cutting plane mode flag.
     *
     * @param cuttingPlaneMode
     */
    public void setCuttingPlaneMode(boolean cuttingPlaneMode) {
        this.cuttingPlaneMode = cuttingPlaneMode;
    }

    public boolean isCuttingPlaneMode() {
        return cuttingPlaneMode;
    }

    /**
     * Sets shading mode flag.
     *
     * @param shadingMode
     */
    public void setShadingMode(boolean shadingMode) {
        this.shadingMode = shadingMode;
    }

    public RaycastRendererPanel getPanel() {
        return panelFront;
    }

    public TransferFunction2DEditor getTF2DPanel() {
        return tfEditor2DFront;
    }

    public TransferFunctionEditor getTFPanel() {
        return tfEditor;
    }

    public TransferFunction2DEditor getTF2DPanelBack() {
        return tfEditor2DBack;
    }

    public TransferFunctionEditor getTFPanelBack() {
        return tfEditorBack;
    }

    //////////////////////////////////////////////////////////////////////
    /////////////////// PRIVATE FUNCTIONS AND ATTRIBUTES /////////////////
    //////////////////////////////////////////////////////////////////////
    /**
     * OpenGL View Matrix. The shape (4x4) remains constant.
     */
    private final double[] _viewMatrix = new double[4 * 4];

    /**
     * Vector used to draw the cutting plane.
     */
    private final double[] _planeU = new double[3];

    /**
     * Vector used to draw the cutting plane.
     */
    private final double[] _planeV = new double[3];

    /**
     * Do NOT modify.
     *
     * Draws the bounding box around the volume.
     *
     * @param gl OpenGL API.
     */
    private void drawBoundingBox(GL2 gl) {
        gl.glPushAttrib(GL2.GL_CURRENT_BIT);
        gl.glDisable(GL2.GL_LIGHTING);
        gl.glColor4d(1.0, 1.0, 1.0, 1.0);
        gl.glLineWidth(1.5f);
        gl.glEnable(GL.GL_LINE_SMOOTH);
        gl.glHint(GL.GL_LINE_SMOOTH_HINT, GL.GL_NICEST);
        gl.glEnable(GL.GL_BLEND);
        gl.glBlendFunc(GL.GL_SRC_ALPHA, GL.GL_ONE_MINUS_SRC_ALPHA);

        gl.glBegin(GL.GL_LINE_LOOP);
        gl.glVertex3d(-volume.getDimX() / 2.0, -volume.getDimY() / 2.0, volume.getDimZ() / 2.0);
        gl.glVertex3d(-volume.getDimX() / 2.0, volume.getDimY() / 2.0, volume.getDimZ() / 2.0);
        gl.glVertex3d(volume.getDimX() / 2.0, volume.getDimY() / 2.0, volume.getDimZ() / 2.0);
        gl.glVertex3d(volume.getDimX() / 2.0, -volume.getDimY() / 2.0, volume.getDimZ() / 2.0);
        gl.glEnd();

        gl.glBegin(GL.GL_LINE_LOOP);
        gl.glVertex3d(-volume.getDimX() / 2.0, -volume.getDimY() / 2.0, -volume.getDimZ() / 2.0);
        gl.glVertex3d(-volume.getDimX() / 2.0, volume.getDimY() / 2.0, -volume.getDimZ() / 2.0);
        gl.glVertex3d(volume.getDimX() / 2.0, volume.getDimY() / 2.0, -volume.getDimZ() / 2.0);
        gl.glVertex3d(volume.getDimX() / 2.0, -volume.getDimY() / 2.0, -volume.getDimZ() / 2.0);
        gl.glEnd();

        gl.glBegin(GL.GL_LINE_LOOP);
        gl.glVertex3d(volume.getDimX() / 2.0, -volume.getDimY() / 2.0, -volume.getDimZ() / 2.0);
        gl.glVertex3d(volume.getDimX() / 2.0, -volume.getDimY() / 2.0, volume.getDimZ() / 2.0);
        gl.glVertex3d(volume.getDimX() / 2.0, volume.getDimY() / 2.0, volume.getDimZ() / 2.0);
        gl.glVertex3d(volume.getDimX() / 2.0, volume.getDimY() / 2.0, -volume.getDimZ() / 2.0);
        gl.glEnd();

        gl.glBegin(GL.GL_LINE_LOOP);
        gl.glVertex3d(-volume.getDimX() / 2.0, -volume.getDimY() / 2.0, -volume.getDimZ() / 2.0);
        gl.glVertex3d(-volume.getDimX() / 2.0, -volume.getDimY() / 2.0, volume.getDimZ() / 2.0);
        gl.glVertex3d(-volume.getDimX() / 2.0, volume.getDimY() / 2.0, volume.getDimZ() / 2.0);
        gl.glVertex3d(-volume.getDimX() / 2.0, volume.getDimY() / 2.0, -volume.getDimZ() / 2.0);
        gl.glEnd();

        gl.glBegin(GL.GL_LINE_LOOP);
        gl.glVertex3d(-volume.getDimX() / 2.0, volume.getDimY() / 2.0, -volume.getDimZ() / 2.0);
        gl.glVertex3d(-volume.getDimX() / 2.0, volume.getDimY() / 2.0, volume.getDimZ() / 2.0);
        gl.glVertex3d(volume.getDimX() / 2.0, volume.getDimY() / 2.0, volume.getDimZ() / 2.0);
        gl.glVertex3d(volume.getDimX() / 2.0, volume.getDimY() / 2.0, -volume.getDimZ() / 2.0);
        gl.glEnd();

        gl.glBegin(GL.GL_LINE_LOOP);
        gl.glVertex3d(-volume.getDimX() / 2.0, -volume.getDimY() / 2.0, -volume.getDimZ() / 2.0);
        gl.glVertex3d(-volume.getDimX() / 2.0, -volume.getDimY() / 2.0, volume.getDimZ() / 2.0);
        gl.glVertex3d(volume.getDimX() / 2.0, -volume.getDimY() / 2.0, volume.getDimZ() / 2.0);
        gl.glVertex3d(volume.getDimX() / 2.0, -volume.getDimY() / 2.0, -volume.getDimZ() / 2.0);
        gl.glEnd();

        gl.glDisable(GL.GL_LINE_SMOOTH);
        gl.glDisable(GL.GL_BLEND);
        gl.glEnable(GL2.GL_LIGHTING);
        gl.glPopAttrib();
    }

    /**
     * Do NOT modify.
     *
     * Draws the cutting plane through.
     *
     * @param gl OpenGL API.
     */
    private void drawCuttingPlane(GL2 gl) {
        gl.glPushAttrib(GL2.GL_CURRENT_BIT);
        gl.glDisable(GL2.GL_LIGHTING);
        gl.glColor4d(1.0, 1.0, 1.0, 1.0);
        gl.glLineWidth(2f);
        gl.glEnable(GL.GL_LINE_SMOOTH);
        gl.glHint(GL.GL_LINE_SMOOTH_HINT, GL.GL_NICEST);
        gl.glEnable(GL.GL_BLEND);
        gl.glBlendFunc(GL.GL_SRC_ALPHA, GL.GL_ONE_MINUS_SRC_ALPHA);

        double D = Math.sqrt(Math.pow(volume.getDimX(), 2) + Math.pow(volume.getDimY(), 2) + Math.pow(volume.getDimZ(), 2)) / 2;

        gl.glBegin(GL.GL_LINE_LOOP);
        gl.glVertex3d(-_planeU[0] * D - _planeV[0] * D, -_planeU[1] * D - _planeV[1] * D, -_planeU[2] * D - _planeV[2] * D);
        gl.glVertex3d(_planeU[0] * D - _planeV[0] * D, _planeU[1] * D - _planeV[1] * D, _planeU[2] * D - _planeV[2] * D);
        gl.glVertex3d(_planeU[0] * D + _planeV[0] * D, _planeU[1] * D + _planeV[1] * D, _planeU[2] * D + _planeV[2] * D);
        gl.glVertex3d(-_planeU[0] * D + _planeV[0] * D, -_planeU[1] * D + _planeV[1] * D, -_planeU[2] * D + _planeV[2] * D);
        gl.glEnd();

        gl.glDisable(GL.GL_LINE_SMOOTH);
        gl.glDisable(GL.GL_BLEND);
        gl.glEnable(GL2.GL_LIGHTING);
        gl.glPopAttrib();
    }

    /**
     * Do NOT modify this function.
     *
     * Sets the Iso value.
     *
     * @param newColor
     */
    public void setIsoValueFront(float isoValueFront) {
        this.isoValueFront = isoValueFront;
    }

    /**
     * Do NOT modify this function.
     *
     * Sets the Iso value.
     *
     * @param newColor
     */
    public void setIsoValueBack(float isoValueBack) {
        this.isoValueBack = isoValueBack;
    }

    /**
     * Do NOT modify this function.
     *
     * Sets the Iso Color.
     *
     * @param newColor
     */
    public void setIsoColorFront(TFColor newColor) {
        this.isoColorFront.r = newColor.r;
        this.isoColorFront.g = newColor.g;
        this.isoColorFront.b = newColor.b;
    }

    /**
     * Do NOT modify this function.
     *
     * Sets the Iso Color.
     *
     * @param newColor
     */
    public void setIsoColorBack(TFColor newColor) {
        this.isoColorBack.r = newColor.r;
        this.isoColorBack.g = newColor.g;
        this.isoColorBack.b = newColor.b;
    }

    /**
     * Do NOT modify this function.
     *
     * Resets the image with 0 values.
     */
    private void resetImage() {
        for (int j = 0; j < image.getHeight(); j++) {
            for (int i = 0; i < image.getWidth(); i++) {
                image.setRGB(i, j, 0);
            }
        }
    }

    /**
     * Do NOT modify this function.
     *
     * Computes the increments according to sample step and stores the result in
     * increments.
     *
     * @param increments Vector to store the result.
     * @param rayVector Ray vector.
     * @param sampleStep Sample step.
     */
    private void computeIncrementsB2F(double[] increments, double[] rayVector, double sampleStep) {
        // we compute a back to front compositing so we start increments in the oposite direction than the pixel ray
        VectorMath.setVector(increments, -rayVector[0] * sampleStep, -rayVector[1] * sampleStep, -rayVector[2] * sampleStep);
    }

    /**
     * Do NOT modify this function.
     *
     * Packs a color into a Integer.
     *
     * @param r Red component of the color.
     * @param g Green component of the color.
     * @param b Blue component of the color.
     * @param a Alpha component of the color.
     * @return
     */
    private static int computePackedPixelColor(double r, double g, double b, double a) {
        int c_alpha = a <= 1.0 ? (int) Math.floor(a * 255) : 255;
        int c_red = r <= 1.0 ? (int) Math.floor(r * 255) : 255;
        int c_green = g <= 1.0 ? (int) Math.floor(g * 255) : 255;
        int c_blue = b <= 1.0 ? (int) Math.floor(b * 255) : 255;
        int pixelColor = (c_alpha << 24) | (c_red << 16) | (c_green << 8) | c_blue;
        return pixelColor;
    }

    /**
     * Do NOT modify this function.
     *
     * Computes the entry and exit of a view vector with respect the faces of
     * the volume.
     *
     * @param p Point of the ray.
     * @param viewVec Direction of the ray.
     * @param entryPoint Vector to store entry point.
     * @param exitPoint Vector to store exit point.
     */
    private void computeEntryAndExit(double[] p, double[] viewVec, double[] entryPoint, double[] exitPoint) {

        for (int i = 0; i < 3; i++) {
            entryPoint[i] = -1;
            exitPoint[i] = -1;
        }

        double[] plane_pos = new double[3];
        double[] plane_normal = new double[3];
        double[] intersection = new double[3];

        VectorMath.setVector(plane_pos, volume.getDimX(), 0, 0);
        VectorMath.setVector(plane_normal, 1, 0, 0);
        intersectFace(plane_pos, plane_normal, p, viewVec, intersection, entryPoint, exitPoint);

        VectorMath.setVector(plane_pos, 0, 0, 0);
        VectorMath.setVector(plane_normal, -1, 0, 0);
        intersectFace(plane_pos, plane_normal, p, viewVec, intersection, entryPoint, exitPoint);

        VectorMath.setVector(plane_pos, 0, volume.getDimY(), 0);
        VectorMath.setVector(plane_normal, 0, 1, 0);
        intersectFace(plane_pos, plane_normal, p, viewVec, intersection, entryPoint, exitPoint);

        VectorMath.setVector(plane_pos, 0, 0, 0);
        VectorMath.setVector(plane_normal, 0, -1, 0);
        intersectFace(plane_pos, plane_normal, p, viewVec, intersection, entryPoint, exitPoint);

        VectorMath.setVector(plane_pos, 0, 0, volume.getDimZ());
        VectorMath.setVector(plane_normal, 0, 0, 1);
        intersectFace(plane_pos, plane_normal, p, viewVec, intersection, entryPoint, exitPoint);

        VectorMath.setVector(plane_pos, 0, 0, 0);
        VectorMath.setVector(plane_normal, 0, 0, -1);
        intersectFace(plane_pos, plane_normal, p, viewVec, intersection, entryPoint, exitPoint);
    }

    /**
     * Do NOT modify this function.
     *
     * Checks if a line intersects a plane.
     *
     * @param plane_pos Position of plane.
     * @param plane_normal Normal of plane.
     * @param line_pos Position of line.
     * @param line_dir Direction of line.
     * @param intersection Vector to store intersection.
     * @return True if intersection happens. False otherwise.
     */
    private static boolean intersectLinePlane(double[] plane_pos, double[] plane_normal,
            double[] line_pos, double[] line_dir, double[] intersection) {

        double[] tmp = new double[3];

        for (int i = 0; i < 3; i++) {
            tmp[i] = plane_pos[i] - line_pos[i];
        }

        double denom = VectorMath.dotproduct(line_dir, plane_normal);
        if (Math.abs(denom) < 1.0e-8) {
            return false;
        }

        double t = VectorMath.dotproduct(tmp, plane_normal) / denom;

        for (int i = 0; i < 3; i++) {
            intersection[i] = line_pos[i] + t * line_dir[i];
        }

        return true;
    }

    /**
     * Do NOT modify this function.
     *
     * Checks if it is a valid intersection.
     *
     * @param intersection Vector with the intersection point.
     * @param xb
     * @param xe
     * @param yb
     * @param ye
     * @param zb
     * @param ze
     * @return
     */
    private static boolean validIntersection(double[] intersection, double xb, double xe, double yb,
            double ye, double zb, double ze) {

        return (((xb - 0.5) <= intersection[0]) && (intersection[0] <= (xe + 0.5))
                && ((yb - 0.5) <= intersection[1]) && (intersection[1] <= (ye + 0.5))
                && ((zb - 0.5) <= intersection[2]) && (intersection[2] <= (ze + 0.5)));

    }

    /**
     * Do NOT modify this function.
     *
     * Checks the intersection of a line with a plane and returns entry and exit
     * points in case intersection happens.
     *
     * @param plane_pos Position of plane.
     * @param plane_normal Normal vector of plane.
     * @param line_pos Position of line.
     * @param line_dir Direction of line.
     * @param intersection Vector to store the intersection point.
     * @param entryPoint Vector to store the entry point.
     * @param exitPoint Vector to store the exit point.
     */
    private void intersectFace(double[] plane_pos, double[] plane_normal,
            double[] line_pos, double[] line_dir, double[] intersection,
            double[] entryPoint, double[] exitPoint) {

        boolean intersect = intersectLinePlane(plane_pos, plane_normal, line_pos, line_dir,
                intersection);
        if (intersect) {

            double xpos0 = 0;
            double xpos1 = volume.getDimX();
            double ypos0 = 0;
            double ypos1 = volume.getDimY();
            double zpos0 = 0;
            double zpos1 = volume.getDimZ();

            if (validIntersection(intersection, xpos0, xpos1, ypos0, ypos1,
                    zpos0, zpos1)) {
                if (VectorMath.dotproduct(line_dir, plane_normal) < 0) {
                    entryPoint[0] = intersection[0];
                    entryPoint[1] = intersection[1];
                    entryPoint[2] = intersection[2];
                } else {
                    exitPoint[0] = intersection[0];
                    exitPoint[1] = intersection[1];
                    exitPoint[2] = intersection[2];
                }
            }
        }
    }

    /**
     * Do NOT modify this function.
     *
     * Calculates the pixel coordinate for the given parameters.
     *
     * @param pixelCoord Vector to store the result.
     * @param volumeCenter Location of the center of the volume.
     * @param uVec uVector.
     * @param vVec vVector.
     * @param i Pixel i.
     * @param j Pixel j.
     */
    private void computePixelCoordinatesFloat(double pixelCoord[], double volumeCenter[], double uVec[], double vVec[], float i, float j) {
        // Coordinates of a plane centered at the center of the volume (volumeCenter and oriented according to the plane defined by uVec and vVec
        float imageCenter = image.getWidth() / 2;
        pixelCoord[0] = uVec[0] * (i - imageCenter) + vVec[0] * (j - imageCenter) + volumeCenter[0];
        pixelCoord[1] = uVec[1] * (i - imageCenter) + vVec[1] * (j - imageCenter) + volumeCenter[1];
        pixelCoord[2] = uVec[2] * (i - imageCenter) + vVec[2] * (j - imageCenter) + volumeCenter[2];
    }

    /**
     * Do NOT modify this function.
     *
     * Same as
     * {@link RaycastRenderer#computePixelCoordinatesFloat(double[], double[], double[], double[], float, float)}
     * but for integer pixel coordinates.
     *
     * @param pixelCoord Vector to store the result.
     * @param volumeCenter Location of the center of the volume.
     * @param uVec uVector.
     * @param vVec vVector.
     * @param i Pixel i.
     * @param j Pixel j.
     */
    private void computePixelCoordinates(double pixelCoord[], double volumeCenter[], double uVec[], double vVec[], int i, int j) {
        // Coordinates of a plane centered at the center of the volume (volumeCenter and oriented according to the plane defined by uVec and vVec
        int imageCenter = image.getWidth() / 2;
        pixelCoord[0] = uVec[0] * (i - imageCenter) + vVec[0] * (j - imageCenter) + volumeCenter[0];
        pixelCoord[1] = uVec[1] * (i - imageCenter) + vVec[1] * (j - imageCenter) + volumeCenter[1];
        pixelCoord[2] = uVec[2] * (i - imageCenter) + vVec[2] * (j - imageCenter) + volumeCenter[2];
    }

    /**
     * Do NOT modify this function.
     *
     * Calculates the pixel coordinate for the given parameters. It calculates
     * the coordinate having the center (0,0) of the view plane aligned with the
     * center of the volume and moved a distance equivalent to the diagonal to
     * make sure we are far enough.
     *
     * @param pixelCoord Vector to store the result.
     * @param viewVec View vector (ray).
     * @param uVec uVector.
     * @param vVec vVector.
     * @param i Pixel i.
     * @param j Pixel j.
     */
    private void computePixelCoordinatesBehindFloat(double pixelCoord[], double viewVec[], double uVec[], double vVec[], float i, float j) {
        int imageCenter = image.getWidth() / 2;
        // Pixel coordinate is calculate having the center (0,0) of the view plane aligned with the center of the volume and moved a distance equivalent
        // to the diaganal to make sure I am far away enough.

        double diagonal = Math.sqrt((volume.getDimX() * volume.getDimX()) + (volume.getDimY() * volume.getDimY()) + (volume.getDimZ() * volume.getDimZ())) / 2;
        pixelCoord[0] = uVec[0] * (i - imageCenter) + vVec[0] * (j - imageCenter) + viewVec[0] * diagonal + volume.getDimX() / 2.0;
        pixelCoord[1] = uVec[1] * (i - imageCenter) + vVec[1] * (j - imageCenter) + viewVec[1] * diagonal + volume.getDimY() / 2.0;
        pixelCoord[2] = uVec[2] * (i - imageCenter) + vVec[2] * (j - imageCenter) + viewVec[2] * diagonal + volume.getDimZ() / 2.0;
    }

    /**
     * Do NOT modify this function.
     *
     * Same as
     * {@link RaycastRenderer#computePixelCoordinatesBehindFloat(double[], double[], double[], double[], int, int)}
     * but for integer pixel coordinates.
     *
     * @param pixelCoord Vector to store the result.
     * @param viewVec View vector (ray).
     * @param uVec uVector.
     * @param vVec vVector.
     * @param i Pixel i.
     * @param j Pixel j.
     */
    private void computePixelCoordinatesBehind(double pixelCoord[], double viewVec[], double uVec[], double vVec[], int i, int j) {
        int imageCenter = image.getWidth() / 2;
        // Pixel coordinate is calculate having the center (0,0) of the view plane aligned with the center of the volume and moved a distance equivalent
        // to the diaganal to make sure I am far away enough.

        double diagonal = Math.sqrt((volume.getDimX() * volume.getDimX()) + (volume.getDimY() * volume.getDimY()) + (volume.getDimZ() * volume.getDimZ())) / 2;
        pixelCoord[0] = uVec[0] * (i - imageCenter) + vVec[0] * (j - imageCenter) + viewVec[0] * diagonal + volume.getDimX() / 2.0;
        pixelCoord[1] = uVec[1] * (i - imageCenter) + vVec[1] * (j - imageCenter) + viewVec[1] * diagonal + volume.getDimY() / 2.0;
        pixelCoord[2] = uVec[2] * (i - imageCenter) + vVec[2] * (j - imageCenter) + viewVec[2] * diagonal + volume.getDimZ() / 2.0;
    }
}
