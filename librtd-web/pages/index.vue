<template>
  <div>
    <textarea v-model="inputSequence" cols="30" rows="10"></textarea>
    <input type="number" v-model="kval" min="1">
    {{ inputSequence ? result : '' }}
  </div>
</template>

<script lang="ts">
import Vue from 'vue'

export default Vue.extend({
  data() {
    return {
      inputSequence: '',
      kval: 2,
    }
  },
  computed: {
    result() {
      //try {
        return sameKmerReturnTimeDistribution(this.inputSequence, this.kval)
      // } catch (error) {
      //   return error.toString()
      // }
    },
  },
})
</script>
