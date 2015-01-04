import os
import shutil
import gzip
import urllib
from urlparse import urljoin

EVENTS_ROOT = 'http://www.psi.utoronto.ca/~hannes/bento-seq-events/'

AS_EVENTS = {
    'hg19': urljoin(EVENTS_ROOT, 'hg19.refseq.event_set.bento-seq.tab.gz'),
    'hg38': urljoin(EVENTS_ROOT, 'hg38.event_set.bento-seq.tab.gz'),
    'mm9': urljoin(EVENTS_ROOT, 'mm9.refseq.event_set.bento-seq.tab.gz'),
    'mm10': urljoin(EVENTS_ROOT, 'mm10.refseq.event_set.bento-seq.tab.gz')
}

DATA_HOME = os.path.abspath(
    os.path.expanduser(
        os.environ.get(
            'BENTOSEQ_HOME',
            os.path.join('~', '.bento-seq'))))

def get_data_home():
    if not os.path.isdir(DATA_HOME):
        os.makedirs(DATA_HOME)
    # XXX: ensure it is dir and readable
    return DATA_HOME


def set_data_home(newpath):
    global DATA_HOME
    DATA_HOME = newpath
    return get_data_home()


def clear_data_home():
    """Delete all the content of the data home cache."""
    data_home = get_data_home()
    shutil.rmtree(data_home)

def count_lines(filename, comment_char=None):
    if filename.endswith('gz') or filename.endswith('gzip'):
        f = gzip.open(filename)
    else:
        f = open(filename)

    if comment_char is None:
        lines = sum(1 for _ in f)
    else:
        return sum(1 for line in f if not line.startswith(comment_char))

    return lines

def fetch(genome):
    if not os.path.isdir(get_data_home()):
        os.makedirs(get_data_home())

    try:
        url = AS_EVENTS[genome]
    except NameError:
        raise NameError("Unknown genome: %s" % genome)

    dest = os.path.join(get_data_home(), os.path.basename(url))
    try:
        gzip.open(dest, 'rb').close()
    except IOError:
        print "Downloading %s from %s." % (genome, url)
        downloader = urllib.urlopen(url)
        data = downloader.read()
        tmp = open(dest, 'wb')
        tmp.write(data)
        tmp.close()
        gzip.open(dest, 'rb').close()

    return dest

def open_event_file(event_string):
    if event_string in AS_EVENTS:
        event_string = fetch(event_string)

    if event_string.endswith('gz') or event_string.endswith('gzip'):
        f = gzip.open(event_string, 'rb')
    else:
        f = open(event_string, 'rb')

    return f