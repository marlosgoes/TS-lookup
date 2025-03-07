function hout = dm_boxplot_1(x,percentiles,colouri)
  
  if (size(percentiles,1)~=1)
    error('percentiles must be a row vector')
  end

  if (size(percentiles,2)~=7)
    error('# columns in percentiles must be 7')
  end
  
  if length(colouri) == 1;
      for ii=1:size(x,2);
    colour{ii} = colouri;
      end
  else
        colour = colouri;
  end
  
  y1 = prctile(x,percentiles(1));
  y2 = prctile(x,percentiles(2));
  y3 = prctile(x,percentiles(3));
  y4 = prctile(x,percentiles(4));
  y5 = prctile(x,percentiles(5));
  y6 = prctile(x,percentiles(6));
  y7 = prctile(x,percentiles(7));
  y_mean = nanmean(x);
  
  box_width = 1/2;%size(x,2)/18;
  whisker_width = box_width/2;
  
  for i=1:size(x,2);
    x_centre(i) = i;
    clear dm_colour;
    
  switch colour{i}
      case 'y'
   dm_colour =  [1 1 0];
	
      case 'm'
dm_colour = [1 0 1] ;

      case 'c'
dm_colour = [0 1 1] ;
     
      case 'r'

dm_colour = [1 0 0] ;
   
      case 'g'
dm_colour = [0 1 0] ;
      case 'b'
dm_colour = [0 0 1] ;
     
      case 'w'

dm_colour = [1 1 1] ;

      case 'k'
dm_colour = [0 0 0] ;

      case 'o'
      dm_colour = [1 0.7 0];    
      
      otherwise 
      dm_colour = colour(i);    
  end



%     if strcmp(colour(i),'o')% colour(i) == 'o'
%       dm_colour = [1 0.7 0];
%     else
%       dm_colour = colour(i);
%     end
    plot(x_centre(i),y1(i),'ow','markerfacecolor',dm_colour,'markersize',5)
    hold on;
    plot([x_centre(i)-whisker_width/2,x_centre(i)+whisker_width/2],[ ... 
        y2(i),y2(i)],'-','color',dm_colour) 
    plot([x_centre(i),x_centre(i)],[y2(i),y3(i)],'-','color',dm_colour)
      
    oi= ciplot([y3(i),y3(i)],[y5(i),y5(i)],...
         [x_centre(i)-box_width/2,x_centre(i)+box_width/2],...
          dm_colour);set(oi,'facealpha',.3)
      
     plot([x_centre(i)-box_width/2,x_centre(i)+box_width/2,x_centre(i)+box_width/2, ...
          x_centre(i)-box_width/2,x_centre(i)-box_width/2], ...
         [y3(i),y3(i),y5(i),y5(i),y3(i)],'-','color',dm_colour) 
  
     
    plot([x_centre(i)-box_width/2,x_centre(i)+box_width/2],[y4(i),y4(i)], ...
         'color',dm_colour) 
    plot([x_centre(i),x_centre(i)],[y5(i),y6(i)],'color',dm_colour) 
    plot([x_centre(i)-whisker_width/2,x_centre(i)+whisker_width/2],[y6(i),y6(i)], ...
         'color',dm_colour) 
    %plot(x_centre(i),y7(i),'color',dm_colour) 
    plot(x_centre(i),y7(i),'ow','markerfacecolor',dm_colour,'markersize',5) 
    plot(x_centre(i),y_mean(i),'*','color',dm_colour,'markersize',5) 
  end
  
  set(gca,'XTick',x_centre)
  set(gca,'XLim',[0.5 0.5+size(x,2)])
  
  return
  
  
