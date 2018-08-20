# Use an official Python runtime as a parent image
FROM ubuntu:14.04
FROM python:2.7

# Set the working directory to /app
WORKDIR /

# Copy the current directory contents into the container at /app
ADD . /

# Install any needed packages specified in requirements.txt
ENV PATH="$PATH:/root/.local/bin"
RUN pip install --user numpy==1.11.2
RUN pip install --user scipy==0.13.3
RUN pip install --user matplotlib
RUN pip install --user astropy==1.0.3
RUN pip install --user cosmolopy
RUN pip install --user pyfits==3.3
RUN pip install --user Pillow
#RUN pip install --user PIL  --allow-unverified PIL --allow-all-external
#RUN pip install --trusted-host pypi.python.org -r requirements.txt

# Make port 80 available to the world outside this container
EXPOSE 8000

# Define environment variable
#ENV NAME World

# Run app.py when the container launches
#CMD ["./galfit"]
RUN chmod +x /mocks.py
ENTRYPOINT ["/mocks.py"]
CMD [""]
#CMD ["python", "mocks.py"]

