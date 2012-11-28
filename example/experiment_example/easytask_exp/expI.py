#!/usr/bin/env python
# -*- encoding: utf-8 -*-
# expI.py
#
# (c) 2012 Dominik Wabersich <dominik.wabersich [aet] gmail.com>
# GPL 3.0+ or (cc) by-sa (http://creativecommons.org/licenses/by-sa/3.0/)
#
# created 2012-10-30
# last mod 2012-10-31 11:15 DW
#


## CONSTANTS ######
KEY1 = 'a'
KEY2 = 'l'
ISI = 0.2 # inter stimulus interval
BSI = 0.2 # before stimulus interval
FI = 0.2 # feedback intervall
NSTIM = 100
SCREENSIZE = [1280,1024]
################### 

from psychopy import visual, misc, core, gui, event
from random import shuffle
import time, os

class Exp:

    def __init__(self, win, fixation):
        ''' Initialize Exp object
            * bsi -- wait before stimulus interval
            * isi -- inter stimulus interval
            * fi -- feedback interval
        '''
        self.win = win
        self.fixation = fixation
        self.key1 = KEY1
        self.key2 = KEY2
        self.isi = ISI
        self.bsi = BSI
        self.fi = FI

    def show_trial(self, vis_stim, give_fb=False):
        ''' show one trial
            * vis_stim - stimulus to show
        '''
        response = False
        self.fixation.draw()
        self.win.flip()
        core.wait(self.bsi)
        vis_stim.draw()
        timekeeper = core.Clock() # start stopwatch
        event.clearEvents()
        self.win.flip()
        while True:
            keys = event.getKeys(timeStamped=timekeeper)
            for key in keys:
                if (key[0] == self.key2 and (vis_stim.name == 'dark' or vis_stim.name == 'undef')):
                    feedback = visual.TextStim(self.win, text='CORRECT', color=[-1,1,-1], height=1)
                    response = (key[0], str(1), str(key[1])[:8])
                elif (key[0] == self.key1 and (vis_stim.name == 'bright' or vis_stim.name == 'undef')):
                    feedback = visual.TextStim(self.win, text='CORRECT', color=[-1,1,-1], height=1)
                    response = (key[0], str(1), str(key[1])[:8])
                elif (key[0] == self.key2):
                    feedback = visual.TextStim(self.win, text='WRONG', color=[1,-1,-1], height=1)
                    response = (key[0], str(0), str(key[1])[:8])
                elif (key[0] == self.key1):
                    feedback = visual.TextStim(self.win, text='WRONG', color=[1,-1,-1], height=1)
                    response = (key[0], str(0), str(key[1])[:8])
                elif (key[0] == 'escape'):
                    core.quit()
            if response:
                break
        if (give_fb):
            feedback.draw()
            self.win.flip()
            core.wait(self.fi)
        self.win.flip()
        core.wait(self.isi)
        return(response)

    def show_inst(self, text='Press any key to continue...'):
        ''' function to show text on the screen for instructions
            * text - text to show on screen
        '''
        instruction = visual.TextStim(self.win, text=text, color=[1,1,1],
            height=1)
        response = False
        instruction.draw()
        event.clearEvents()
        self.win.flip()
        while True:
            keys = event.getKeys()
            for key in keys:
                if (key == self.key1 or key == self.key2):
                    response = True
                elif (key == 'escape'):
                    core.quit()
            if response:
                break
        self.win.flip()
        core.wait(self.isi)


if __name__ == '__main__':
    # show dialog box for patricipant info
    try:
        exp_info = misc.fromFile('last_expinfo.pickle')
    except:
        exp_info = {'Nr':'01', 'Age (18-99)':22, 'Sex (m/f)':'m', 
                'Hand (r/l)':'r'}
    exp_info_gui = gui.DlgFromDict(exp_info)
    if exp_info_gui.OK:
        misc.toFile('last_expinfo.pickle', exp_info)
    else:
        core.quit()
    exp_info['Timestamp'] = time.strftime('%Y-%m-%d-%H%M%S')

    # create window 
    win = visual.Window(SCREENSIZE, monitor='labscreen', units='cm',
            color=[-1,-1,-1], allowGUI=False)
    # create fixation cross
    fixation = visual.ShapeStim(win=win, pos=[0,0], 
                vertices=((-.02,.02),(-.02,.2),(.02,.2),(.02,.02),
                    (.2,.02),(.2,-.02),(.02,-.02),(.02,-.2),(-.02,-.2),
                    (-.02,-.02),(-.2,-.02),(-.2,.02)),
                lineColor=None,fillColor=(1,1,1), name='fix')
    # create stim list
    vis_stims1 = []
    for i in range(int(NSTIM/3)):
        vis_stims1.append(visual.PatchStim(win=win, mask='circle', pos=[0,0],
            color=[0,0,0], name='dark', tex=None, size=10))
    for i in range(int(NSTIM/3)):
        vis_stims1.append(visual.PatchStim(win=win, mask='circle', pos=[0,0],
            color=[0.5,0.5,0.5], name='undef', tex=None, size=10))
    for i in range(int(NSTIM/3)):
        vis_stims1.append(visual.PatchStim(win=win, mask='circle', pos=[0,0],
            color=[1,1,1], name='bright', tex=None, size=10))

    # create exp instance
    exp = Exp(win, fixation)

    # check if data directory does exist
    try: 
        os.stat('data/')
    except OSError:
        os.mkdir('data') # make data directory if it does not exist
    try: 
        os.stat('data/readme_exp1.txt')
    except OSError:
        with open('data/readme_exp1.txt', 'w') as infofile:
            infofile.write('Experiment 1: Easytask' 
                    + '--------------------------------------\n'
                    + '\nInfo about data goes here')

############################################################################
    # draw a bit...
    exp.fixation.draw()
    exp.win.flip()
    stim_text = visual.TextStim(exp.win, text='AAA', color=[1,1,1], height=1)
    stim_text.draw()
    exp.win.flip()
    patch_stim = visual.PatchStim(exp.win, mask='circle', pos=[0,0],
        color=[0,0,0], tex=None, size=10)
    patch_stim.draw()
    exp.win.flip()
    core.wait(0.1)
    exp.win.flip()
############################################################################


    # open file for datawriting
    with open('data/exp1_' + str(exp_info['Nr']) + '_' 
            + str(exp_info['Sex (m/f)'])
            + str(exp_info['Age (18-99)'])
            + str(exp_info['Hand (r/l)']) + '_'
            + str(exp_info['Timestamp'])
            + '.txt', 'w') as datafile:
        datafile.write('Participant: ' + str(exp_info['Nr']) + '\n')
        datafile.write('Sex: ' + str(exp_info['Sex (m/f)']) + '\n')
        datafile.write('Age: ' + str(exp_info['Age (18-99)']) + '\n')
        datafile.write('Hand: ' + str(exp_info['Hand (r/l)']) + '\n')
        datafile.write('Timestamp: ' + str(exp_info['Timestamp']) + '\n')
        datafile.write('\n------------------------ Data -------------------------\n\n')
        datafile.write('%2.2s%7.4s%7.3s%3.1s%10.2s\n' 
                %('Id', 'Stim', 'Key', 'C', 'RT'))
    
        # instructions:
        exp.show_inst('Easy Task:\n'
                      + 'There is one bright '
                      + 'and one dark circle. '
                      + 'There is also a mixture circle '
                      + 'between bright and dark.\n'
                      + 'Press the left reaction key, '
                      + 'if the circle belongs to the bright class, '
                      + 'and the right reaction key if the circle belongs '
                      + 'to the dark class. The mixed circle belongs to '
                      + 'both classes, so press either the left or the '
                      + 'right reaction key.\n'
                      + '\nPress any key to continue...')

        # training block, 20 trials
        exp.show_inst('Easy Task.\n'
                      + '\nPress any key to start training block')

        shuffle(vis_stims1)
        for vis_stim in vis_stims1[:20]:
            participant_response = exp.show_trial(vis_stim, True)
            datafile.write('%2.2s' %exp_info['Nr'].encode('utf8', 'ignore'))
            datafile.write('%7.5s' %(vis_stim.name)) # Stim
            datafile.write('%7.5s%3.1s%10.8s\n' 
                %participant_response)

        # exp block 1
        exp.show_inst('Training done.\n'
                      + '\nPress any key to start experimental block 1/2')

        shuffle(vis_stims1)
        for vis_stim in vis_stims1:
            participant_response = exp.show_trial(vis_stim)
            datafile.write('%2.2s' %exp_info['Nr'].encode('utf8', 'ignore'))
            datafile.write('%7.5s' %(vis_stim.name)) # Stim
            datafile.write('%7.5s%3.1s%10.8s\n' 
                %participant_response)

        # exp block 2
        exp.show_inst('Block 1 done.\n'
                      + '\nPress any key to start experimental block 2/2')

        shuffle(vis_stims1)
        for vis_stim in vis_stims1:
            participant_response = exp.show_trial(vis_stim)
            datafile.write('%2.2s' %exp_info['Nr'].encode('utf8', 'ignore'))
            datafile.write('%7.5s' %(vis_stim.name)) # Stim
            datafile.write('%7.5s%3.1s%10.8s\n' 
                %participant_response)
        
        # done
        exp.show_inst('Done!')
    
    exp.win.close()
    core.quit()
