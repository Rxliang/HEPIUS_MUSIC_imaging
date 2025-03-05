function c = pinmap_ge408

% from verasonics GE connector pdf

c=[];

c.ch = [];
c.ind = [];
i=1;
c.ind{i} = 159:-1:128;
c.ch{i}.letter = 'A';
c.ch{i}.num = 2:33;

i=i+1;
c.ind{i} = 190:-2:160;
c.ch{i}.letter = 'B';
c.ch{i}.num = 3:2:33;

i=i+1;
c.ind{i} = [222:-2:192];
c.ch{i}.letter = 'C';
c.ch{i}.num = 2:2:32;

i=i+1;
c.ind{i} = [191:-2:161];
c.ch{i}.letter = 'C';
c.ch{i}.num = 3:2:33;

i=i+1;
c.ind{i} = [223:-2:193];
c.ch{i}.letter = 'D';
c.ch{i}.num = 3:2:33;

i=i+1;
c.ind{i} = [255:-1:224];
c.ch{i}.letter = 'E';
c.ch{i}.num = 2:1:33;

i=i+1;
c.ind{i} = [31:-1:0];
c.ch{i}.letter = 'F';
c.ch{i}.num = 2:1:33;

i=i+1;
c.ind{i} = [62:-2:32];
c.ch{i}.letter = 'G';
c.ch{i}.num = 2:2:32;

i=i+1;
c.ind{i} = [94:-2:64];
c.ch{i}.letter = 'H';
c.ch{i}.num = 2:2:32;

i=i+1;
c.ind{i} = [63:-2:33];
c.ch{i}.letter = 'H';
c.ch{i}.num = 3:2:33;

i=i+1;
c.ind{i} = [95:-2:65];
c.ch{i}.letter = 'J';
c.ch{i}.num = 2:2:32;

i=i+1;
c.ind{i} = [127:-1:96];
c.ch{i}.letter = 'K';
c.ch{i}.num = 2:1:33;


c.pad = [];
c.padChannel = [];
c.mapping = [];
c.channelNo = [];

cnt=0;

for i = 1:length(c.ch)
  numEl = length(c.ind{i});
  for j = 1:numEl
    chan=1+c.ind{i}(j);
    c.pad{chan} = [c.ch{i}.letter num2str(c.ch{i}.num(j))];
    c.padChannel(chan) = c.ch{i}.num(j);    
    cnt=cnt+1;
    %mapping{cnt,1} = [pad{ind{i}(j)} ' ' num2str(Trans.HVMux.Aperture(cnt))];
    %ch6s(cnt) = Trans.HVMux.Aperture(cnt);
    c.mapping{cnt,1} = [c.pad{chan} ' ' num2str(chan)];
    c.channelNo(cnt)=chan;
  end
 
end

% from computeTrans: mapping to Verasonics logical channels: these
% numbers are logical numbers, the index are physical channels on
% the GE 

c.Connector1D = [ 97    98    99   100   101   102   103   104   105   106   111   108   109   110   107   112, ...
                    171   180   170   179   169   178   168   177   164   176   163   175   162   174   161   173, ...
                    120   121   116   122   115   123   114   124   113   125   119   126   118   127   117   128, ...
                    181   192   182   191   183   190   184   189   165   188   166   187   167   186   172   185, ...
                    33    44    34    43    35    42    36    41    37    48    38    47    39    46    40    45, ...
                    241   245   242   246   243   247   244   248   228   229   227   230   226   231   225   232, ...
                    49    57    51    58    53    59    55    61    50    63    52    60    54    62    56    64, ...
                    236   240   235   239   234   238   233   237   252   256   251   255   250   254   249   253, ...
                    8    12     7    11     6    10     5     9     4    16     3    15     2    14     1    13, ...
                    205   197   206   198   193   199   194   200   195   209   196   210   207   211   208   212, ...
                    17    28    18    27    19    26    20    25    21    32    22    31    23    30    24    29, ...
                    216   224   215   223   214   222   213   221   204   217   203   218   202   219   201   220, ...
                    83    84    82    85    81    86    80    87    68    69    67    70    66    71    65    72, ...
                    145   152   146   151   147   150   148   149   129   136   130   135   131   134   132   133, ...
                    92    96    91    95    90    94    89    93    88    79    75    78    74    77    73    76, ...
                    153   157   154   158   155   159   156   160   137   141   138   142   139   143   140   144 ]';
                
