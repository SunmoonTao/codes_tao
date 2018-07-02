
class DNA:
    def __init__(self, sequence, DNAtype='Linear'):
        self.seq=sequence
        if DNAtype=='Circular':
            self.cyclic=True
        else:
            self.cyclic=False
        self.sequence_length=len(sequence)
        if self.cyclic:
            self.type='Circular'
        else:
            self.type='Linear'
        
    def getdna(self,start, stop):
        if start > 0 \
            and stop > 0 \
            and start <= stop \
            and stop <= self.sequence_length:
                return self.seq[start:stop]
        if self.cyclic:
            if start > self.sequence_length:
                start = start-self.sequence_length
            if stop > self.sequence_length:
                stop = stop-self.sequence_length
            if start <= 0:
                start = self.sequence_length + start
            if stop <= 0:
                stop = self.sequence_length + stop
            if start > stop:
                return self.seq[start:] + self.seq[:stop]
            return self.seq[start:stop]
        raise ValueError('indices out of bounds')       