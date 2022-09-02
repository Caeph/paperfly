import sys
import openstack

inputfile = sys.argv[1]
storedname = sys.argv[2]
container = sys.argv[3]

conn = openstack.connect(cloud='openstack')
obj = conn.object_store.create_object(container=container,
                                      name=storedname,
                                      data=open(inputfile, 'r'),
                                      content_type='txt')