from abc import abstractmethod, ABC

import numpy as np
from skimage import draw, filters
from matplotlib import pyplot as plt
import pydicom
from pydicom.dataset import Dataset, FileMetaDataset


def clip_add(image1, image2, dtype=np.uint16):
    """Clip the image to the dtype extrema. Otherwise the bits will flip."""
    return np.clip(image1 + image2, np.iinfo(dtype).min, np.iinfo(dtype).max).astype(dtype)


def even_round(num):
    num = int(round(num))
    return num + num % 2


class Layer(ABC):

    @abstractmethod
    def apply(self, image, pixel_size):
        pass


class ConeLayer(Layer):

    def __init__(self, cone_size_mm=10, cax_offset_mm=(0, 0), alpha=1.0):
        self.cone_size_mm = cone_size_mm
        self.cax_offset_mm = cax_offset_mm
        self.alpha = alpha

    def apply(self, image: np.ndarray, pixel_size):
        cone_size_pix = (self.cone_size_mm / 2) / pixel_size
        cax_offset_pix = [x / pixel_size + (shape/2 - 0.5) for x, shape in zip(self.cax_offset_mm, image.shape)]
        rr, cc = draw.disk(cax_offset_pix, cone_size_pix, shape=image.shape)
        rr = np.round(rr).astype(np.int)
        cc = np.round(cc).astype(np.int)
        image[rr, cc] = int(np.iinfo(image.dtype).max * self.alpha)
        return image


class FieldLayer(Layer):

    def __init__(self, field_size_mm=(10, 10), cax_offset_mm=(0, 0), alpha=1.0):
        self.field_size_mm = field_size_mm
        self.cax_offset_mm = cax_offset_mm
        self.alpha = alpha

    def apply(self, image: np.ndarray, pixel_size):
        field_size_pix = [even_round(f / pixel_size) for f in self.field_size_mm]
        field_start = [x / pixel_size + (shape/2) - field_size/2 for x, shape, field_size in zip(self.cax_offset_mm, image.shape, field_size_pix)]
        field_end = [x / pixel_size + (shape/2) + field_size/2 - 1 for x, shape, field_size in zip(self.cax_offset_mm, image.shape, field_size_pix)]
        # -1 due to skimage implementation of [start:(end+1)]
        rr, cc = draw.rectangle(field_start, end=field_end, shape=image.shape)
        rr = np.round(rr).astype(np.int)
        cc = np.round(cc).astype(np.int)
        image[rr, cc] = int(np.iinfo(image.dtype).max * self.alpha)
        return image


class BBLayer(ConeLayer):

    def __init__(self, bb_size_mm=5, cax_offset_mm=(0, 0), alpha=0.5):
        super().__init__(cone_size_mm=bb_size_mm, cax_offset_mm=cax_offset_mm, alpha=alpha)


class GaussianLayer(Layer):

    def __init__(self, sigma_mm=2):
        self.sigma_mm = sigma_mm

    def apply(self, image, pixel_size):
        sigma_pix = self.sigma_mm / pixel_size
        return filters.gaussian(image, sigma_pix, preserve_range=True).astype(image.dtype)


class RandomNoiseLayer(Layer):

    def __init__(self, mean=0.0, sigma=0.01):
        self.mean = mean
        self.sigma = sigma

    def apply(self, image, pixel_size):
        noise = np.random.normal(self.mean, self.sigma, size=image.shape)
        return clip_add(image, noise, dtype=image.dtype)


class ConstantLayer(Layer):

    def __init__(self, constant):
        self.constant = constant

    def apply(self, image, pixel_size):
        return clip_add(image, self.constant, dtype=image.dtype)


class AS500Image:
    pixel_size = 0.784
    shape = (384, 512)

    def __init__(self):
        self.image = np.zeros(self.shape, np.uint16)

    def add_layer(self, layer):
        self.image = layer.apply(self.image, self.pixel_size)

    def generate_dicom(self, file_out_name: str, gantry_angle: float = 0.0, coll_angle: float = 0.0, table_angle: float = 0.0, sid: float = 1500.0):
        file_meta = FileMetaDataset()
        # Main data elements
        ds = Dataset()
        ds.ImageType = ['DERIVED', 'SECONDARY', 'PORTAL']
        ds.SOPClassUID = '1.2.840.10008.5.1.4.1.1.481.1'
        ds.SOPInstanceUID = '1.2.840.113854.141883099300381770008774160544352783139.1.1'
        ds.StudyDate = '20150212'
        ds.ContentDate = '20150212'
        ds.StudyTime = '124120'
        ds.ContentTime = '124120'
        ds.AccessionNumber = ''
        ds.Modality = 'RTIMAGE'
        ds.ConversionType = 'WSD'
        ds.Manufacturer = 'IMPAC Medical Systems, Inc.'
        ds.ReferringPhysicianName = ''
        ds.StudyDescription = 'QA'
        ds.SeriesDescription = 'D + Gantry_0'
        ds.PhysiciansOfRecord = 'Awesome Physician'
        ds.OperatorsName = ''
        ds.ManufacturerModelName = 'MOSAIQ'
        ds.PatientName = 'Lutz^Test Tool'
        ds.PatientID = 'zzzBAC_Lutz'
        ds.PatientBirthDate = ''
        ds.SoftwareVersions = '2.41.01J0'
        ds.StudyInstanceUID = '1.2.840.113854.141883099300381770008774160544352783139'
        ds.SeriesInstanceUID = '1.2.840.113854.141883099300381770008774160544352783139.1'
        ds.StudyID = '348469'
        ds.SeriesNumber = "4597199"
        ds.InstanceNumber = "0"
        ds.PatientOrientation = ''
        ds.PositionReferenceIndicator = ''
        ds.SamplesPerPixel = 1
        ds.PhotometricInterpretation = 'MONOCHROME2'
        ds.Rows = self.image.shape[0]
        ds.Columns = self.image.shape[1]
        ds.BitsAllocated = 16
        ds.BitsStored = 16
        ds.HighBit = 15
        ds.PixelRepresentation = 0
        ds.RTImageLabel = 'D'
        ds.RTImagePlane = 'NORMAL'
        ds.XRayImageReceptorAngle = "0.0"
        ds.ImagePlanePixelSpacing = [self.pixel_size, self.pixel_size]
        ds.RTImagePosition = [-200.70400, 150.52800]
        ds.RadiationMachineSAD = "1000.0"
        ds.RTImageSID = str(sid)
        ds.PrimaryDosimeterUnit = 'MU'
        ds.GantryAngle = str(gantry_angle)
        ds.BeamLimitingDeviceAngle = str(coll_angle)
        ds.PatientSupportAngle = str(table_angle)
        ds.PixelData = self.image  # XXX Array of 393216 bytes excluded

        ds.file_meta = file_meta
        ds.is_implicit_VR = True
        ds.is_little_endian = True
        ds.save_as(file_out_name, write_like_original=False)


class AS1000Image(AS500Image):
    pixel_size = 0.392
    shape = (768, 1024)

    def generate_dicom(self, file_out_name: str, gantry_angle: float = 0.0, coll_angle: float = 0.0, table_angle: float = 0.0, sid: float = 1500.0):
        # File meta info data elements
        file_meta = FileMetaDataset()

        # Main data elements
        ds = Dataset()
        ds.ImageType = ['DERIVED', 'SECONDARY', 'PORTAL']
        ds.SOPClassUID = '1.2.840.10008.5.1.4.1.1.481.1'
        ds.SOPInstanceUID = '1.2.840.113854.323870129946883845737820671794195198978.1.1'
        ds.StudyDate = '20140819'
        ds.ContentDate = '20140819'
        ds.StudyTime = '130944'
        ds.ContentTime = '130944'
        ds.AccessionNumber = ''
        ds.Modality = 'RTIMAGE'
        ds.ConversionType = 'WSD'
        ds.Manufacturer = 'IMPAC Medical Systems, Inc.'
        ds.ReferringPhysicianName = ''
        ds.StudyDescription = 'QA'
        ds.SeriesDescription = 'Q + Couch_270'
        ds.PhysiciansOfRecord = 'I am king'
        ds.OperatorsName = ''
        ds.ManufacturerModelName = 'MOSAIQ'
        ds.PatientName = 'Albert Einstein'
        ds.PatientID = 'abc123'
        ds.PatientBirthDate = ''
        ds.SoftwareVersions = '2.41.01J0'
        ds.StudyInstanceUID = '1.2.840.113854.323870129946883845737820671794195198978'
        ds.SeriesInstanceUID = '1.2.840.113854.323870129946883845737820671794195198978.1'
        ds.StudyID = '348469'
        ds.SeriesNumber = "4290463"
        ds.InstanceNumber = "0"
        ds.PatientOrientation = ''
        ds.PositionReferenceIndicator = ''
        ds.SamplesPerPixel = 1
        ds.PhotometricInterpretation = 'MONOCHROME2'
        ds.Rows = self.shape[0]
        ds.Columns = self.shape[1]
        ds.BitsAllocated = 16
        ds.BitsStored = 16
        ds.HighBit = 15
        ds.PixelRepresentation = 0
        ds.RTImageLabel = 'Q'
        ds.RTImagePlane = 'NORMAL'
        ds.XRayImageReceptorAngle = "0.0"
        ds.ImagePlanePixelSpacing = [self.pixel_size, self.pixel_size]
        ds.RTImagePosition = [-200.70400, 150.523400]
        ds.RadiationMachineSAD = "1000.0"
        ds.RTImageSID = str(sid)
        ds.PrimaryDosimeterUnit = 'MU'
        ds.GantryAngle = str(gantry_angle)
        ds.BeamLimitingDeviceAngle = str(coll_angle)
        ds.PatientSupportAngle = str(table_angle)
        # ds.CurveDimensions = 2
        # ds.NumberOfPoints = 4
        # ds.TypeOfData = 'ROI'
        # ds.CurveDescription = 'VContour 30'
        # ds.AxisUnits = ['PIXL', 'PIXL']
        # ds.DataValueRepresentation = 3
        # ds.CurveLabel = 'Field Edge (Q:MV)'
        # ds.CurveData = b'\x00\x00\x00\x00 .\x81@\x00\x00\x00\x00\x10\\z@\x00\x00\x00\x00 .\x81@\x00\x00\x00\x00\x90\x93u@\x00\x00\x00\x00\xc0\x93}@\x00\x00\x00\x00\x90\x93u@\x00\x00\x00\x00\xc0\x93}@\x00\x00\x00\x00\x10\\z@'
        ds.PixelData = self.image # XXX Array of 1572864 bytes excluded

        ds.file_meta = file_meta
        ds.is_implicit_VR = True
        ds.is_little_endian = True
        ds.save_as(file_out_name, write_like_original=False)


class AS1200Image(AS500Image):
    pixel_size = 0.336
    shape = (1280, 1280)

    def generate_dicom(self, file_out_name: str, gantry_angle: float = 0.0, coll_angle: float = 0.0, table_angle: float = 0.0, sid: float = 1500.0):
        file_meta = FileMetaDataset()
        file_meta.FileMetaInformationGroupLength = 196
        file_meta.FileMetaInformationVersion = b'\x00\x01'
        file_meta.MediaStorageSOPClassUID = '1.2.840.10008.5.1.4.1.1.481.1'
        file_meta.MediaStorageSOPInstanceUID = '1.2.246.352.64.1.5468686515961995030.4457606667843517571'
        file_meta.TransferSyntaxUID = '1.2.840.10008.1.2'
        file_meta.ImplementationClassUID = '1.2.246.352.70.2.1.120.1'
        file_meta.ImplementationVersionName = 'MergeCOM3_410'

        # Main data elements
        ds = Dataset()
        ds.SpecificCharacterSet = 'ISO_IR 100'
        ds.ImageType = ['ORIGINAL', 'PRIMARY', 'PORTAL']
        ds.InstanceCreationDate = '20161230'
        ds.InstanceCreationTime = '215510'
        ds.SOPClassUID = '1.2.840.10008.5.1.4.1.1.481.1'
        ds.SOPInstanceUID = '1.2.246.352.64.1.5468686515961995030.4457606667843517571'
        ds.StudyDate = '20161230'
        ds.ContentDate = '20161230'
        ds.StudyTime = '215441.936'
        ds.ContentTime = '215441.936'
        ds.AccessionNumber = ''
        ds.Modality = 'RTIMAGE'
        ds.ConversionType = ''
        ds.Manufacturer = 'Varian Medical Systems'
        ds.ReferringPhysicianName = ''
        ds.StationName = 'NDS-WKS-SN9999'
        ds.OperatorsName = 'King Kong'
        ds.ManufacturerModelName = 'VMS.XI.Service'
        ds.PatientName = 'Grace Hopper'
        ds.PatientID = 'VMS.XI.Service'
        ds.PatientBirthDate = '19000101'
        ds.PatientBirthTime = '000000'
        ds.PatientSex = ''
        ds.SoftwareVersions = '2.5.13.2'
        ds.StudyInstanceUID = '1.2.246.352.64.4.5644626269434644263.1905029945372990626'
        ds.SeriesInstanceUID = '1.2.246.352.64.2.5508761605912087323.11665958260371685307'
        ds.StudyID = 'fdd794f2-8520-4c4a-aecc-e4446c1730ff'
        ds.SeriesNumber = None
        ds.AcquisitionNumber = "739774555"
        ds.InstanceNumber = "1"
        ds.PatientOrientation = ''
        ds.FrameOfReferenceUID = '1.2.246.352.64.3.4714322356925391886.9391210174715030407'
        ds.PositionReferenceIndicator = ''
        ds.SamplesPerPixel = 1
        ds.PhotometricInterpretation = 'MONOCHROME1'
        ds.PlanarConfiguration = 0
        ds.Rows = self.shape[0]
        ds.Columns = self.shape[1]
        ds.BitsAllocated = 16
        ds.BitsStored = 16
        ds.HighBit = 15
        ds.PixelRepresentation = 0
        ds.WindowCenter = "32767.0"
        ds.WindowWidth = "65535.0"
        ds.RescaleIntercept = "0.0"
        ds.RescaleSlope = "1.0"
        ds.RescaleType = 'US'
        ds.RTImageLabel = 'MV_180'
        ds.RTImageDescription = ""
        ds.ReportedValuesOrigin = 'ACTUAL'
        ds.RTImagePlane = 'NORMAL'
        ds.XRayImageReceptorTranslation = [0.00, 0.00, 1000-sid]
        ds.XRayImageReceptorAngle = "0.0"
        ds.RTImageOrientation = [1, 0, 0, 0, -1, 0]
        ds.ImagePlanePixelSpacing = [self.pixel_size, self.pixel_size]
        ds.RTImagePosition = [-214.872, 214.872]
        ds.RadiationMachineName = 'TrueBeam from Hell'
        ds.RadiationMachineSAD = "1000.0"
        ds.RTImageSID = str(sid)

        ds.PrimaryDosimeterUnit = 'MU'
        ds.GantryAngle = str(gantry_angle)
        ds.BeamLimitingDeviceAngle = str(coll_angle)
        ds.PatientSupportAngle = str(table_angle)
        ds.TableTopVerticalPosition = "-24.59382842824"
        ds.TableTopLongitudinalPosition = "200.813502948597"
        ds.TableTopLateralPosition = "3.00246706215532"
        ds.PixelData = self.image  # XXX Array of 3276800 bytes excluded

        ds.file_meta = file_meta
        ds.is_implicit_VR = True
        ds.is_little_endian = True
        ds.save_as(file_out_name, write_like_original=False)
