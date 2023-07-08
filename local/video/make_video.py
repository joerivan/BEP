import cv2
import os
from PIL import Image, ImageOps
import numpy as np

def main():

    # set path of this folder
    folder = os.path.join(os.path.dirname(__file__))
    
    # clean up images and stich together
    for img in np.arange(1,301):
        plot_image = Image.open(os.path.join(folder, 'plot_images', '%i.png' %img)).convert("RGBA")
        
        #cleaning up system image
        syst_image = Image.open(os.path.join(folder, 'syst_images', '%i.png' %img)).convert("RGBA")
        line_width = 2
        border = (line_width,line_width,line_width,line_width)
        border_img = ImageOps.expand(syst_image, border=border, fill='#000000')
        scaled_img = border_img.resize((2000,2000))
        
        image = stich_images(plot_image, scaled_img)
        resize_image = image.resize((4096,2731))
        resize_image.save(os.path.join(folder, 'images', '%i.png' %img))#, quality=90)
    
    
    # make video
    video_name = os.path.join(os.path.dirname(__file__),'ideo_name.mp4')
    
    images = []
    for file in np.arange(1,301):
        images.append(f'{file}.png')
        
    frame = cv2.imread(os.path.join(folder,'images', images[0]))
    
    height, width, layers = frame.shape
    fourcc = cv2.VideoWriter_fourcc(*'avc1')
    
    video = cv2.VideoWriter(video_name, fourcc, 30, (width,height))

    # adding images to video
    for image in images:
        video.write(cv2.imread(os.path.join(folder, 'images', image)))
    
    # change order images and add again to video
    images.reverse()
    
    for image in images:
        video.write(cv2.imread(os.path.join(folder, 'images', image)))
    
    # finish video
    cv2.destroyAllWindows()
    video.release()
    
def stich_images(plot_img,syst_img):
    """
    stiches two images together
    """
    plot_img.paste(syst_img, (4200,1850), mask=syst_img)
    
    return plot_img
    
    
if __name__ == '__main__':
    main()
