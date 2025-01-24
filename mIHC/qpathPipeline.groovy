// set colors (unmixed images)
setChannelNames(
    'DAPI',
    'Ki67',
    'CD45',
    'PDL2',
    'CK',
    'PDL1',
    'Tissue'
)

// create pixel annotation
setImageType('FLUORESCENCE');
resetSelection();
createAnnotationsFromPixelClassifier("isTissue", 0.0, 0.0)

// set whole image as ROI
createFullImageAnnotation(true);

// stardist with cell expansion
import qupath.ext.stardist.StarDist2D

// Specify the model file (you will need to change this!)
def pathModel = '/Users/Dru/Documents/Andreas/Bordet/Projects/rejuveron/qtiffs/quPathCellSegModels/dsb2018_heavy_augment.pb'

def stardist = StarDist2D.builder(pathModel)
        .threshold(0.5)              // Probability (detection) threshold
        .channels('DAPI')            // Select detection channel
        .normalizePercentiles(1, 99) // Percentile normalization
        .pixelSize(0.5)              // Resolution for detection
        .cellExpansion(10.0)         // Approximate cells based upon nucleus expansion
        .measureShape()              // Add shape measurements
        .measureIntensity()          // Add cell measurements (in all compartments)
        .includeProbability(true)    // Add probability as a measurement (enables later filtering)
        .build()
        
// Run detection for the selected objects
def imageData = getCurrentImageData()
def pathObjects = getSelectedObjects()
if (pathObjects.isEmpty()) {
    Dialogs.showErrorMessage("StarDist", "Please select a parent object!")
    return
}
stardist.detectObjects(imageData, pathObjects)
println 'Done!'

// classify
runObjectClassifier("all_classes")

// get image path
name = GeneralTools.stripExtension(getCurrentImageData().getServer().getMetadata().getName())
// save measurements
saveDetectionMeasurements(buildFilePath("/Users/Dru/Downloads/", name+".txt"))
