# R code for the study
# Pollination interactions reveal direct costs and indirect benefits of plant–plant facilitation for ecosystem engineers
# Gianalberto Losapio & Christian Schöb
# Journal of Plant Ecology
# October 2019

library(reshape)
library(nlme)
library(car)
library(lsmeans)
library(lmerTest)
library(effects)

#
insect.record_cb <- read.csv("insect_recordR.csv", head=T)
insect.record_cb <- insect.record_cb[which(insect.record_cb$pl.code=="horspi"|insect.record_cb$pl.code=="aretet"),]

insect.record_cb <- insect.record_cb[which(insect.record_cb$in.guild=="pollinator"),]

seeds_cb <- read.csv("/seeds_cb.csv", head=T)

# 1. n flowers

df_fl <- melt(sndata_cb, id.vars=c("microh","plot","block","plsp"), measure.vars = c("aretet.fl","horspi.fl"))
df_fl <- df_fl[-c(29:84),]
rownames(df_fl) <- NULL


### 2. n. visitors/flower

sndata_cb$insects <- 0
sndata_cb$poll <- 0
sndata_cb$nonpol <- 0

for(i in 1:nrow(sndata_cb)){
	sndata_cb$insects[i] <- nrow(insect.record_cb[which(insect.record_cb$foundation==sndata_cb$fs[i]&insect.record_cb$microh==sndata_cb$microh[i]&insect.record_cb$plot==sndata_cb$plot[i]),])
	sndata_cb$poll[i] <- nrow(insect.record_cb[which(insect.record_cb$foundation==sndata_cb$fs[i]&insect.record_cb$microh==sndata_cb$microh[i]&insect.record_cb$plot==sndata_cb$plot[i]&insect.record_cb$in.guild=="pollinator"),])
	sndata_cb$nonpol[i] <- sndata_cb$insects[i]-sndata_cb$poll[i]
}

head(sndata_cb)

df_costbenef <- melt(sndata_cb, measure.vars = c("aretet.fl","horspi.fl"))
df_costbenef <- df_costbenef[-c(29:84),]
colnames(df_costbenef)[12] <- "fsflowers"
df_costbenef$variable <- NULL

df_costbenef$pollfl <- df_costbenef$poll/df_costbenef$fsflowers

### add diversity of pollinators
insect.record_cb$polsp <- paste(insect.record_cb$genus,insect.record_cb$species,sep=".")
df_costbenef$polrich<-0

for(i in 1:nrow(df_costbenef)){
df_costbenef$polrich[i]<-length(unique(insect.record_cb$polsp[which(insect.record_cb$foundation==df_costbenef$fs[i]&insect.record_cb$microh==df_costbenef$microh[i]&insect.record_cb$plot==df_costbenef$plot[i])]))
}

df_costbenef$polrichfl <- df_costbenef$polrich/df_costbenef$fsflowers

## reproductive success

df_costbenef <- merge(df_costbenef,seeds_cb, by.x=c("fs","microh","plot"), by.y=c("foundation","microh","plot"))

df_costbenef$fruflo <- df_costbenef$fruits/df_costbenef$fsflowers

df_costbenef$seedsfr <- df_costbenef$seeds/df_costbenef$fruits

df_costbenef$seedsfr <- df_costbenef$seeds/df_costbenef$fruits

df_costbenef$seedsflo <- df_costbenef$seeds/df_costbenef$fsflowers

# 1 number of flowers (per plant)
mod.fsfl <-glmer(fsflowers ~microh*fs + (1 |block), data= df_costbenef, family="poisson")

Anova(mod.fsfl, test="Chisq", type=3)
summary(mod.fsfl)

summary(Effect("microh", mod.fsfl))
round((236.3009-385.4235)/236.3009*100)

summary(allEffects(mod.fsfl))

round((432.3190-206.7347)/432.3190*100)
round((343.6150-270.0956)/343.6150*100)

plot(allEffects(mod.fsfl))

# effect size
fixef(mod.fsfl)[-1]/(summary(mod.fsfl)$coefficients[-1,2]*sqrt(56))

# contrast
emmeans(mod.fsfl, pairwise ~ microh*fs)


# 2 pollinator visitation rate
mod.pollfl <-lmer(pollfl ~microh*fs + (1 |block), data= df_costbenef)

Anova(mod.pollfl, test="Chisq")
summary(mod.pollfl)

summary(Effect("microh", mod.pollfl))

round((0.01906148-0.01222546)/0.01906148*100)

plot(allEffects(mod.pollfl))

fixef(mod.pollfl)[-1]/(summary(mod.pollfl)$coefficients[-1,2]*sqrt(56))

# 2b pollinator richness
mod.pollsp <-lmer(polrichfl ~microh*fs + (1 |block), data= df_costbenef)

Anova(mod.pollsp, test="Chisq")
summary(mod.pollsp)

summary(Effect("microh", mod.pollsp))

round((0.015037031-0.008536726)/0.015037031*100)

plot(allEffects(mod.pollsp))

emmeans(mod.pollsp, pairwise ~ microh)

fixef(mod.pollsp)[-1]/(summary(mod.pollsp)$coefficients[-1,2]*sqrt(56))

## 3 reproductive success

### fruit set
cor.test(df_costbenef$pollfl, df_costbenef$polrichfl)

mod.fruflo.null<-lmer(fruflo ~ 1+ (1 |block), data= df_costbenef)

mod.fruflo.1<- lmer(fruflo ~ microh*fs+ (1 |block), data= df_costbenef)

mod.fruflo.2a<-lmer(fruflo ~ microh*fs+ pollfl + (1 |block), data= df_costbenef)

mod.fruflo.2b<-lmer(fruflo ~ microh*fs+ polrichfl + (1 |block), data= df_costbenef)

mod.fruflo.2c<-lmer(fruflo ~ microh*fs+ poll + (1 |block), data= df_costbenef)

anova(mod.fruflo.null,mod.fruflo.1,mod.fruflo.2a,mod.fruflo.2b,mod.fruflo.2c)
anova(mod.fruflo.null,mod.fruflo.1,mod.fruflo.2b,mod.fruflo.2c,mod.fruflo.2a)
anova(mod.fruflo.null,mod.fruflo.1,mod.fruflo.2c,mod.fruflo.2b,mod.fruflo.2a)

# fruit set final model
mod.fruflo2 <- glmer(cbind(fruits,fl) ~ microh*fs + polrichfl + (1 |block),
                     data= df_costbenef, family=binomial)
Anova(mod.fruflo2, type=3)
summary(mod.fruflo2)

summary(Effect("microh", mod.fruflo2))

summary(Effect("fs", mod.fruflo2))

summary(allEffects(mod.fruflo2))

(0.07485285 + 0.09438946 + 0.11837270 + 0.18221176 + 0.22301444)/5

round((0.07331284-0.11400053)/0.07331284*100)

round((0.11475399-0.09680487)/0.11475399*100)


plot(Effect(c("microh","fs"), mod.fruflo2))

fixef(mod.fruflo2)[-1]/(summary(mod.fruflo2)$coefficients[-1,2]*sqrt(56))

emmeans(mod.fruflo2, pairwise ~ microh* fs)

######
mod.fruits <- glmer(fruits ~ microh*fs+ (1 |block), data= df_costbenef, family=poisson)
Anova(mod.fruits)
summary(mod.fruits)

#### seeds

## 3 a: seeds (per plant)
mod.seeds.a.null <-glmer(seeds ~ 1+ (1 |block), data= df_costbenef, family=poisson)

mod.seeds.a.1 <- glmer(seeds ~ microh*fs+ (1 |block), data= df_costbenef, family=poisson)

mod.seeds.a.2c <- glmer(seeds ~ microh*fs+ poll + (1 |block), data= df_costbenef, family=poisson)

mod.seeds.a.2a <- glmer(seeds ~ microh*fs+ pollfl + (1 |block), data= df_costbenef, family=poisson)

mod.seeds.a.2b <- glmer(seeds ~ microh*fs+ polrichfl + (1 |block), data= df_costbenef, family=poisson)

anova(mod.seeds.a.null,mod.seeds.a.1,mod.seeds.a.2a,mod.seeds.a.2b,mod.seeds.a.2c)
anova(mod.seeds.a.null,mod.seeds.a.1,mod.seeds.a.2b,mod.seeds.a.2c,mod.seeds.a.2a)
anova(mod.seeds.a.null,mod.seeds.a.1,mod.seeds.a.2c,mod.seeds.a.2b,mod.seeds.a.2a)

# seeds per plant final model
mod.seeds<- glmer(seeds ~ microh*fs+ poll + (1 |block), data= df_costbenef, family=poisson)
Anova(mod.seeds, test="Chisq", type=3)
summary(mod.seeds)

summary(Effect("microh", mod.seeds))

summary(allEffects(mod.seeds))

(38.38469 + 42.10309 + 46.18171 +  50.65542 + 55.56252)/5

round((108.33469-65.52511)/108.33469*100)
round((22.74766-15.94080)/22.74766*100)

fixef(mod.seeds)[-1]/(summary(mod.seeds)$coefficients[-1,2]*sqrt(56))

emmeans(mod.seeds, pairwise ~ microh* fs)

## 3 b: seeds per fruit
mod.seeds.b.null <-lmer(seedsfr ~ 1+ (1 |block), data= df_costbenef)

mod.seeds.b.1 <- lmer(seedsfr ~ microh*fs+ (1 |block), data= df_costbenef)

mod.seeds.b.2c <- lmer(seedsfr ~ microh*fs+ poll + (1 |block), data= df_costbenef)

mod.seeds.b.2a <- lmer(seedsfr ~ microh*fs+ pollfl + (1 |block), data= df_costbenef)

mod.seeds.b.2b <- lmer(seedsfr ~ microh*fs+ polrichfl + (1 |block), data= df_costbenef)

anova(mod.seeds.b.null,mod.seeds.b.1,mod.seeds.b.2a,mod.seeds.b.2b,mod.seeds.b.2c)
anova(mod.seeds.b.null,mod.seeds.b.1,mod.seeds.b.2b,mod.seeds.b.2c,mod.seeds.b.2a)
anova(mod.seeds.b.null,mod.seeds.b.1,mod.seeds.b.2c,mod.seeds.b.2b,mod.seeds.b.2a)

# seeds per fruit final model
mod.seedsfru2 <- glmer(cbind(seeds,fruits) ~ microh*fs + (1 |block),
                       data= df_costbenef, family=binomial)
Anova(mod.seedsfru2)
summary(mod.seedsfru2)

## 3 c: seeds per flower
mod.seeds.c.null <-lmer(seedsflo ~ 1+ (1 |block), data= df_costbenef)

mod.seeds.c.1 <- lmer(seedsflo ~ microh*fs+ (1 |block), data= df_costbenef)

mod.seeds.c.2c <- lmer(seedsflo ~ microh*fs+ poll + (1 |block), data= df_costbenef)

mod.seeds.c.2a <- lmer(seedsflo ~ microh*fs+ pollfl + (1 |block), data= df_costbenef)

mod.seeds.c.2b <- lmer(seedsflo ~ microh*fs+ polrichfl + (1 |block), data= df_costbenef)

anova(mod.seeds.c.null,mod.seeds.c.1,mod.seeds.c.2a,mod.seeds.c.2b,mod.seeds.c.2c)
anova(mod.seeds.c.null,mod.seeds.c.1,mod.seeds.c.2b,mod.seeds.c.2c,mod.seeds.c.2a)
anova(mod.seeds.c.null,mod.seeds.c.1,mod.seeds.c.2c,mod.seeds.c.2b,mod.seeds.c.2a)

# seeds per flower final model
mod.seedsflo2 <- glmer(cbind(seeds,fl) ~ microh*fs + pollfl + (1 |block),
                       data= df_costbenef, family=binomial)
Anova(mod.seedsflo2, type=3)
summary(mod.seedsflo2)

summary(allEffects(mod.seedsflo2))

(0.07962603 + 0.14013125 + 0.23487937 + 0.36639293 + 0.44247997)/5

round((0.2640906-0.1683296)/0.2640906*100)
round((22.74766-15.94080)/22.74766*100)

fixef(mod.seedsflo2)[-1]/(summary(mod.seedsflo2)$coefficients[-1,2]*sqrt(56))

emmeans(mod.seedsflo2, pairwise ~ microh* fs)

setEPS
postscript("mod.seedsflo2a.eps", widt=2.5, height=2.5)
plot(Effect(c("microh","fs"),mod.seedsflo2))
dev.off()

pdf("mod.seedsflo2b.pdf", widt=2.5, height=2.5)
plot(Effect(c("pollfl"),mod.seedsflo2))
dev.off()

setEPS()
postscript("modseeds1.eps", widt=2.5, height=2.5)
plot(Effect(c("microh","fs"), mod.seeds))
dev.off()

pdf("modseeds2.pdf", widt=2.5, height=2.5)
plot(Effect(c("poll"), mod.seeds))
dev.off()


###### plots

setEPS()
postscript("fig2a.eps", widt=3.5, height=3.5)
plot(allEffects(mod.fsfl))
dev.off()

setEPS()
postscript("fig2b.eps", widt=3.5, height=3.5)
plot(allEffects(mod.pollfl))
dev.off()

setEPS()
postscript("fig2d.eps", widt=3.5, height=3.5)
plot(allEffects(mod.pollsp))
dev.off()

setEPS()
postscript("fig4c1.eps", widt=2.5, height=2.5)
plot(Effect(c("microh","fs"), mod.fruflo2))
dev.off()

pdf("fig4c3.pdf", widt=2.5, height=2.5)
plot(Effect(c("polrichfl"), mod.fruflo2))
dev.off()
#######

write.csv(df_costbenef, "df_costbenef.csv")

write.csv(round(coef(summary(mod.pollfl)),3), "mod.pollfl.csv")
write.csv(round(coef(summary(mod.fsfl)),3), "mod.fsfl.csv")
write.csv(round(coef(summary(mod.seeds)),3), "mod.seeds.csv")
write.csv(round(coef(summary(mod.seeds.b)),3), "mod.seedsfr.csv")


save.image("losapio_costben_jpe.RData")
