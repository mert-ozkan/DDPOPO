function plot_quest_results(op,int_spdX)

for whPos = 1:8
    for whTrj = 1:2
        for whIll = 1:2
            field_nm = sprintf('Stim%d%d%d',whPos,whTrj,whIll);
            spd(whPos,whTrj,whIll) = op.(field_nm).avg;
        end
    end
end
%%
close all
figure(1);
ax = gca;
for whPos = 1:8
    subplot(1,8,whPos);
    datN = [spd(whPos,1,1), spd(whPos,1,2), spd(whPos,2,1), spd(whPos,1,2)];
    plt = bar([1 2 4 5],datN);
    ylim([min(spd(:))-1,max(spd(:))]+.5);
end

%%
figure(2)
condX = int_spdX(:,1)*100+int_spdX(:,2)*10+int_spdX(:,3);
step = .00001;
x = 0:step:6;
for whPos = 1:8
    for whTrj = 1:2
        for whIll = 1:2
            whSp = (whPos-1)*4+(whTrj-1)*2+whIll;
            
            subplot(8,4,whSp)
            condN = whPos*100+whTrj*10+whIll;
            trl_condN = condX==condN;
            
            x_dat = int_spdX(trl_condN,4);
            y_dat = int_spdX(trl_condN,5);
            x_lim = mean(x_dat) + std(x_dat) * 1 * [-1, 1];
            
            bX = glmfit(x_dat,y_dat,'binomial','link','logit');
            z = bX(1) + (x * bX(2));
            z = 1 ./ (1 + exp(-z));
            plt_dat = scatter(x_dat, y_dat);
            plt_est = line(x,z);
            
            [~,idx] = min(abs(z-.5));
            thresh_logit = x(idx);
            diff((whTrj-1)*2+whIll, whPos) = spd(whPos,whTrj,whIll)-thresh_logit;
            
            if ~isempty(thresh_logit)
                vline(thresh_logit,'r');
            end
            vline(spd(whPos,whTrj,whIll),'g:');
            
            xlim(x_lim);

        end
    end
end
end